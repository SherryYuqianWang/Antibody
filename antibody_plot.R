# Package ######################################################################
## Maths / Biostats -----------------------------------------------------------
library(deSolve) 
library(MASS) #sample from multivariate normal distribution
library(finalfit) #table of regression result
library(survival) #coxph function
## Plotting --------------------------------------------------------------------
library(ggplot2)
library(ggsci)
library(ggh4x)
library(cowplot)
## Data processing -------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(writexl)
library(purrr) #Functional Programming


rm(list=ls())
ABpath <- "~/Documents/UZH/0_random/Sherry/antibody/4. Yun Shan - Antibody Modelling/"
source(paste0(ABpath,"Code/antibody_function.R")) #function for antibody plots

# Setting ######################################################################
Tmin <- 0
Tmax <- 360

step_size <- 1
times <- c(seq(Tmin,Tmax,step_size))

# Population fit ###############################################################
## Import data -----------------------------------------------------------------
pop_fit <- list()
pop_fit[[1]] <- read.csv(paste0(ABpath,"Data/Monolix/r11_wt_igg_03/populationParameters.txt"), row.names = 1)
pop_fit[[2]] <- read.csv(paste0(ABpath,"Data/Monolix/r11_ba1_igg_14/populationParameters.txt"), row.names = 1)
pop_fit[[3]] <- read.csv(paste0(ABpath,"Data/Monolix/r11_wt_iga_07/populationParameters.txt"), row.names = 1)
pop_fit[[4]] <- read.csv(paste0(ABpath,"Data/Monolix/r11_ba1_iga_03/populationParameters.txt"), row.names = 1)

## Generate figure -------------------------------------------------------------
num = 100
ab_list <- c("WT IgG (%)","BA.1 IgG (%)","WT IgA (%)","BA.1 IgA (%)")
plot_settings <- list(
  list(var = "Age", colors = c("#E1BA4E", "#B279B4")),
  list(var = "Gender", colors = c("#D98DAB", "#46C1BE")),
  list(var = "Booster", colors = c("#D55C5C", "#368AC5"))
)

# generate different datasets base on age, gender, and booster type
group_df <- list()

group_age <- create_group_df(1, 1, c(0, 1))
group_gender <- create_group_df(1, c(1, 0), 0)
group_booster <- create_group_df(c(1, 0), 1, 0)

group_df <- list(group_age, group_gender, group_booster)


plot_list = list()

for (i in seq_along(ab_list)) {
  
  pop <- pop_fit[[i]]
  ab <- ab_list[i]
  
  for (j in seq_along(plot_settings)) {
    
    setting <- plot_settings[[j]]
    group_name <- setting$var
    colors <- setting$colors
    
    group_list <- split(group_df[[j]], seq(nrow(group_df[[j]]))) 
    total_AB <- map(group_list, simulate_AB)
    plot_AB <- list_rbind(total_AB) %>%  mutate(Booster = case_when(Booster == 0 ~ "Moderna", 
                                                                    Booster == 1 ~ "Pfizer"), 
                                                Infect = case_when(Infect == 0 ~ "Early Infection",
                                                                   Infect == 1 ~ "Late Infection",
                                                                   Infect == 2 ~ "Uninfected"),
                                                Age = case_when(Age == 0 ~ "<60 years",
                                                                Age == 1 ~ ">=60 years"),
                                                Gender = case_when(Gender == 0 ~ "Female", 
                                                                   Gender == 1 ~ "Male"))
    
    p <- make_plot(data = plot_AB,
                   group_var = group_name,
                   y_label = ab,
                   color_values = colors, 
                   fill_values = colors)
    
    plot_list[[(i - 1) * 3 + j]] <- p
  }
}

  
ggdraw() +
  draw_plot(plot_list[[3]], x = 0.02, y = 0.71, width = 0.44, height = 0.24) +
  draw_plot(plot_list[[6]], x = 0.02, y = 0.50, width = 0.44, height = 0.22) +
  draw_plot(plot_list[[4]], x = 0.02, y = 0.21, width = 0.44, height = 0.24) +
  draw_plot(plot_list[[5]], x = 0.02, y = 0, width = 0.44, height = 0.22) +
  draw_plot(plot_list[[9]], x = 0.47, y = 0.71, width = 0.54, height = 0.24) +
  draw_plot(plot_list[[12]], x = 0.47, y = 0.50, width = 0.54, height = 0.22) +
  draw_plot(plot_list[[10]], x = 0.47, y = 0.21, width = 0.54, height = 0.24) +
  draw_plot(plot_list[[11]], x = 0.47, y = 0, width = 0.54, height = 0.22) +
  draw_plot_label(label = c("A   WT- and BA.1-specific responses:", 
                            "B   BA.1-specific responses:                "), size = 12,
                  x = c(-0.12,-0.12), y = c(0.995, 0.495))

ggsave(paste0(ABpath,"Figures/Figure3_revision1.png"), width = 10.5, height = 7.5,bg = "white")
ggsave(paste0(ABpath,"Figures/Figure3_revision1.pdf"), width = 10.5, height = 7.5,bg = "white")

# Individual fit ###############################################################
## Import data -----------------------------------------------------------------
original <- read.csv(paste0(ABpath,"Data/P360_Sflow_For Keisuke_17Jan2024_remove74.csv"))
original$booster_trans <- strptime(as.character(original$Booster.Vaccination.date), "%d/%m/%Y")
original_cov <- list()
original_cov[[1]] <- original %>% select(Code, booster_trans, DoI, WT.IgG) %>% group_by(Code) %>% slice(2)
original_cov[[2]] <- original %>% select(Code, booster_trans, DoI, BA1.IgG) %>% group_by(Code) %>% slice(2)
original_cov[[3]] <- original %>% select(Code, booster_trans, DoI, WT.IgA) %>% group_by(Code) %>% slice(2)
original_cov[[4]] <- original %>% select(Code, booster_trans, DoI, BA1.IgA) %>% group_by(Code) %>% slice(2)

sg_covid <- list()
for (i in 1:4){
  sg_covid[[i]] <- read.csv(paste0(ABpath,"weekly_cases_per_million.csv"))
  sg_covid[[i]]$date_trans <- strptime(as.character(sg_covid[[i]]$date), "%Y-%m-%d")
}

Estimated <- list()
Estimated[[1]] <- read.csv(paste0(ABpath,"Data/Monolix/r11_wt_igg_03/IndividualParameters/estimatedIndividualParameters.txt"), sep = ",", comment.char = "", header = T)
Estimated[[2]] <- read.csv(paste0(ABpath,"Data/Monolix/r11_ba1_igg_14/IndividualParameters/estimatedIndividualParameters.txt"), sep = ",", comment.char = "", header = T)
Estimated[[3]] <- read.csv(paste0(ABpath,"Data/Monolix/r11_wt_iga_07/IndividualParameters/estimatedIndividualParameters.txt"), sep = ",", comment.char = "", header = T)
Estimated[[4]] <- read.csv(paste0(ABpath,"Data/Monolix/r11_ba1_iga_03/IndividualParameters/estimatedIndividualParameters.txt"), sep = ",", comment.char = "", header = T)

Simulated <- list()
Simulated[[1]] <- read.csv(paste0(ABpath,"Data/Monolix/r11_wt_igg_03/IndividualParameters/simulatedIndividualParameters.txt"), sep = ",", comment.char = "", header = T)
Simulated[[2]] <- read.csv(paste0(ABpath,"Data/Monolix/r11_ba1_igg_14/IndividualParameters/simulatedIndividualParameters.txt"), sep = ",", comment.char = "", header = T)
Simulated[[3]] <- read.csv(paste0(ABpath,"Data/Monolix/r11_wt_iga_07/IndividualParameters/simulatedIndividualParameters.txt"), sep = ",", comment.char = "", header = T)
Simulated[[4]] <- read.csv(paste0(ABpath,"Data/Monolix/r11_ba1_iga_03/IndividualParameters/simulatedIndividualParameters.txt"), sep = ",", comment.char = "", header = T)

## Tables ----------------------------------------------------------------------

#Generate dataset with k1,k2,k3,k4
ind_AB <- list()
ind_AB <- pmap(list(Estimated,original_cov,sg_covid),ind_ab_long)

#Generate dataset with k1,k2
ind_AB <- list()
ind_AB <- map2(Estimated,original_cov,ind_ab_long_cluster)


sheets <- list("wt_igg" = ind_AB[[1]], "ba1_igg" = ind_AB[[2]], "wt_iga" = ind_AB[[3]], "ba1_iga" = ind_AB[[4]])
write_xlsx(sheets, paste0(ABpath,"Output/antibody_individual_protection_80%.xlsx"))

## Figures ---------------------------------------------------------------------
#Generate figures
ind_fit <- list()
ind_fit <- pmap(list(Estimated,Simulated,original_cov),ind_fit_plt)

antibody <-  c("WT IgG","BA1 IgG","WT IgA","BA1 IgA")

for (a in 1:4){
  
  pdf(paste0("plot/",antibody[a], ".pdf"), 11, 10)
  for (i in seq(1, length(unique(ind_fit[[a]]$Code)), 49)) {
    print(
      ggplot(ind_fit[[a]][ind_fit[[a]]$Code %in% levels(ind_fit[[a]]$Code)[i:(i+48)],]) +
        geom_point(aes(x=Day,y=.data[[antibody[a]]]),color="#cc718b",size=1.5, shape=16, stroke = 3) +
        geom_line(aes(x=Day,y=y),lwd=1, color ="#7FA2C5") +
        geom_ribbon(aes(x=Day,ymin=Min,ymax=Max), fill="#7FA2C5", alpha=0.2) +
        facet_wrap(vars(Code), ncol=7, nrow=7)+
        xlab("Day after booster vaccination") +
        ylab(paste0(antibody[a], " binding (%)"))  +
        ylim(0,101)+
        theme(axis.text = element_text(colour = "black"),
              axis.ticks = element_line(colour = "black"),
              axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              legend.position='none',
              axis.title.y = element_text(size=11,family="sans"),
              axis.title.x = element_text(size=11,family="sans")))
    
  }
  dev.off()
}

## Regression ------------------------------------------------------------------
reg_ab <- list()
reg_ab <- map(ind_AB,regression_table)

combine_reg <- map_df(reg_ab, ~as.data.frame(.x), .id = "Antibody") %>%
  mutate(Antibody = case_when(Antibody == "1" ~ "WT IgG",
                              Antibody == "2" ~ "BA1 IgG",
                              Antibody == "3" ~ "WT IgA",
                              Antibody == "4" ~ "BA1 IgA")) %>% select(!(fit_id))

write_xlsx(combine_reg, paste0(ABpath,"Output/regression_result.xlsx"))


# Survival analysis ############################################################
## Data prep -------------------------------------------------------------------
ba1_iga_cat <- read.csv(paste0(ABpath,"ba1_iga_new.csv"))

# prep data frames for survival models
ba1_iga_cat$doi_new[is.na(ba1_iga_cat$doi_new)] <- 339
ba1_iga_df_cat <- filter(ba1_iga_cat, time <= doi_new)
ba1_iga_df_cat$Infect[ba1_iga_df_cat$Infect == 0] <- 1
ba1_iga_df_cat$Infect[ba1_iga_df_cat$Infect == 2] <- 0
#ba1_iga_df_cat <- rename(ba1_iga_df_cat, BA1.IgA = antibody)

time_ind_ba1_iga <- ba1_iga_df_cat %>% select(id, Gender, Infect, age_cat, Booster, day28_ab,
                                              doi_new) %>% distinct()
time_ind_ba1_iga <- tmerge(data1 = time_ind_ba1_iga, data2 = time_ind_ba1_iga, id = id, event = event(doi_new, Infect))
time_dep_ba1_iga <- ba1_iga_df_cat %>% select(id, antibody, time, cases, 
                                              antibody_cat4, day28_cat4)
full_ba1_iga <- tmerge(data1 = time_ind_ba1_iga, data2 = time_dep_ba1_iga, id = id, 
                       antibody = tdc(time, antibody), cases = tdc(time, cases), 
                       cat4 = tdc(time, antibody_cat4), day28_cat4 = tdc(time, day28_cat4))
full_ba1_iga <- filter(full_ba1_iga, tstart != 0)

full_ba1_iga <- full_ba1_iga %>% 
  mutate(case_cat = cases/5920000*100000)

test <- full_ba1_iga %>%
  subset(doi_new==tstop & event==1) %>% arrange(cases)
summary(test$cases)

## Coxph model -----------------------------------------------------------------
ba1_iga_cont <- coxph(Surv(tstart, tstop, event) ~ antibody + 
                        log2(cases) , data = full_ba1_iga)
#+ Gender + age_cat + Booster, data = full_ba1_iga ) 

summary(ba1_iga_cont) # higher ba1 IgA => lower risk of infection


## Prediction ------------------------------------------------------------------
#age_cat <- c(">=60","<60")
#booster <- c("1","0")
#gender <- c("Male","Female")
cases <- c(2653,4347,8099) #weekly cases per million
tstop_values <- seq(10,360,10)
k2_ref <- (exp(log(0.00614)+0.2) + 0.00614)/2
A0_values <- c(1,10,20,30,40,50)

Tmin <- 28
step_size <- 1
tstop_ind <- 388
rep_n = (tstop_ind-Tmin)+1


pars_list <- lapply(A0_values, function(a0) {
  c(A0 = a0, k2 = k2_ref)  # make it a named numeric vector
})



results <- map(pars_list, ~ ode_ab_surv(pars = ., tstop_ind = tstop_ind))
names(results) <- A0_values

results_df <- bind_rows(results, .id = "AO")

pred_df<-data.frame()


for (a in seq_along(A0_values)){
  for (c in 1:3){
  test.dat <- data.frame(antibody = results[[a]][, 2],
                         tstart=c(Tmin:tstop_ind),
                         tstop=c((Tmin+1):(tstop_ind+1)),
                         event=rep(0,rep_n),
                         cases=rep(cases[c],rep_n),
                         A0 = A0_values[a])
                         # age_cat=rep(age_cat[a],rep_n),
                         # Booster=rep(booster[b],rep_n),
  # Gender=rep(gender[g],rep_n))
  pred_df <- rbind(pred_df,test.dat)
  #}
  #}
  }
}

preds <- predict(ba1_iga_cont, newdata = pred_df, type = "survival",se.fit = TRUE)

pred_df$prob <-preds$fit
#pred_df$upr <- preds$fit + (1.96 * preds$se.fit)
#pred_df$lwr <- preds$fit - (1.96 * preds$se.fit)
#pred_df$cases <- as.factor(pred_df$cases)



# Create all combinations
combinations <- expand.grid(A0 = A0_values,
                            cases = cases,
                            tstop = tstop_values)

# Calculate cumulative probability for each combination
pred_df_plt <- combinations %>%
  pmap_dfr(function(A0_val, cases_val, tstop_val) {
     filtered <- pred_df %>%
       filter(A0 == A0_val, cases == cases_val, tstop <= tstop_val)
     
     prob_cum <- prod(filtered$prob)
    
    tibble(A0 = A0_val, cases = cases_val, tstop = tstop_val, prob_cum = prob_cum)
  })


pred_df_plt <- pred_df_plt %>% 
  mutate(#Booster = case_when(Booster == "0" ~ "Moderna",
    #                    Booster == "1" ~ "Pfizer"),
    #age_cat = case_when(age_cat == "<60" ~ "<60 yr",
    #                   age_cat == ">=60" ~ ">=60 yr"),
    cases = case_when(cases == 2653 ~ "Low transmission",
                      cases == 4347 ~ "Intermediate transmission",
                      cases == 8099 ~ "High transmission"),
    cases = factor(cases, levels = c("Low transmission", "Intermediate transmission", "High transmission")),
    tstop_fac = as.factor(tstop),
    A0_fac =as.factor(A0))
    #upr = ifelse(upr<=1,upr,1),
    #lwr = ifelse(lwr<=0,0,lwr))


## Figure ----------------------------------------------------------------------
#ba1_IgA_3 months
low_3m_80 <- 11
mod_6m_80 <- 9
high_3m_80 <- 22

#ba1_IgA_6 months
low_6m_80 <-35
mod_9m_80 <- 25
high_6m_80 <- 42

line_dat1 <- data.frame(
  cases = c("Low transmission", "Low transmission","Low transmission",
            "Intermediate transmission", "Intermediate transmission", "Intermediate transmission",
            "High transmission", "High transmission", "High transmission"),
  x = c(low_3m_80, low_6m_80, 0,
        mod_6m_80, mod_9m_80, 0,
        high_3m_80, high_6m_80, 0),
  xend = c(low_3m_80, low_6m_80, low_6m_80,
           mod_6m_80, mod_9m_80, mod_9m_80,
           high_3m_80, high_6m_80,high_6m_80),
  y = c(0, 0, 80,
        0, 0, 80,
        0, 0, 80),
  yend = c(80, 80, 80,
           80, 80, 80,
           80, 80, 80)
) %>% mutate(
  cases = factor(cases, levels = c("Low transmission", "Intermediate transmission", "High transmission"))
)

figA <- ggplot(data=pred_df_plt %>% filter(tstop %in% c("90","180","270","360")))+
  geom_smooth(aes(x = A0,y = prob_cum*100, linetype=tstop_fac),color="blue",linewidth = 0.5)+
  #geom_ribbon(aes(x = antibody, ymin = lwr*100, ymax = upr*100),alpha=0.5,fill="#003399FF") +
  labs(colour="Group") +
  geom_segment(data = line_dat1, aes(x = x, xend = xend, y = y, yend = yend), color = "darkblue", linetype = "dashed",linewidth=0.3)+
  facet_grid2(~cases,axes = "all")+
  ylim(0,100)+
  xlim(0,51)+
  scale_y_continuous(breaks=seq(0,100,by=20),labels = expression(0,20,40,60,80,100),limits=c(0.0,100)) +
  scale_linetype_manual(values = c("solid","dotted","dotdash","longdash"))+
  xlab("BA.1 IgA binding (%)")+
  ylab("Protection against infection (%)")+
  labs(linetype = "Days after booster vaccination") +  # Change legend title
  theme_classic() +
  theme(legend.position = "bottom")  


#ba1_IgA_3 months
low_20_80 <- 215
mod_20_80 <- 155
high_20_80 <- 95

#ba1_IgA_6 months
low_40_80 <- 290
mod_40_80 <- 225
high_40_80 <- 170

line_dat2 <- data.frame(
  cases = c("Low transmission", "Low transmission","Low transmission",
            "Intermediate transmission", "Intermediate transmission", "Intermediate transmission",
            "High transmission", "High transmission", "High transmission"),
  x = c(low_20_80, low_40_80, 0,
        mod_20_80, mod_40_80, 0,
        high_20_80, high_40_80, 0),
  xend = c(low_20_80, low_40_80, low_40_80,
           mod_20_80, mod_40_80, mod_40_80,
           high_20_80, high_40_80,high_40_80),
  y = c(0, 0, 80,
        0, 0, 80,
        0, 0, 80),
  yend = c(80, 80, 80,
           80, 80, 80,
           80, 80, 80)
) %>% mutate(
  cases = factor(cases, levels = c("Low transmission", "Intermediate transmission", "High transmission"))
)


figB <- ggplot(data=pred_df_plt %>% filter(A0 != "1"))+
  #geom_line(aes(x = tstop,y = prob_cum*100, linetype=A0_fac),color="darkblue")+
  geom_smooth(aes(x = tstop,y = prob_cum*100, linetype=A0_fac),color="blue",linewidth = 0.5)+
  #geom_ribbon(aes(x = antibody, ymin = lwr*100, ymax = upr*100),alpha=0.5,fill="#003399FF") +
  #labs(colour="Group") +
  geom_segment(data = line_dat2, aes(x = x, xend = xend, y = y, yend = yend), color = "darkblue", linetype = "dashed",linewidth=0.3)+
  facet_grid2(~cases,axes = "all")+
  ylim(0,100)+
  xlim(0,360)+
  scale_y_continuous(breaks=seq(0,100,by=20),labels = expression(0,20,40,60,80,100),limits=c(0,100)) +
  scale_linetype_manual(values = c("solid","dotted","dashed","dotdash","longdash"))+
  xlab("Days after booster vaccination")+
  ylab("Protection against infection (%)")+
  labs(linetype = "BA.1 IgA binding (%)") +  # Change legend title
  theme_classic() +
  theme(legend.position = "bottom")  




ggdraw() +
  draw_plot(figA, x = 0.05, y = 0.5, width = 0.9, height = 0.48) +
  draw_plot(figB, x = 0.05, y = 0.0, width = 0.9, height = 0.48) +
  draw_plot_label(label = c("A", 
                            "B"), size = 12,
                  x = c(0.01,0.01), y = c(0.995, 0.495))

ggsave(paste0(ABpath,"Figures/Figure6_revision1.png"), width = 9, height = 6,bg = "white")
ggsave(paste0(ABpath,"Figures/Figure6_revision1.pdf"), width = 9, height = 6,bg = "white")



#Table 2. Summary of binding antibody levels (%)################################
table2_df <- list()
antibody <- c("WT IgG", "BA1 IgG","WT IgA", "BA1 IgA")

for (i in 1:4){
  original_ab <- original[,c(1,10,10+i)]
  colnames(original_ab)[3] <- "ab"
  
  table2_df[[i]] <- original_ab %>% 
    mutate(table2 = case_when(ab < summary(sheets[[i]]$antibody_new)[2] ~ 'low', 
                              ab >= summary(sheets[[i]]$antibody_new)[2] & ab < summary(sheets[[i]]$antibody_new)[5] ~ 'medium',
                              ab >= summary(sheets[[i]]$antibody_new)[5] ~ 'high')) %>% 
    group_by(table2,Day) %>% 
    summarise(mean_IQR = paste0(round(median(ab),digits = 1), " (", 
                                round(quantile(ab, 0.25),digits = 1), ", ", 
                                round(quantile(ab, 0.75),digits = 1), ")"), .groups = "drop") %>%
    mutate(antibody = antibody[i])
}

combine_table2 <- map_df(table2_df, ~as.data.frame(.x)) %>% 
  pivot_wider(names_from = table2,
              values_from = mean_IQR)%>%
  relocate(high, .after = last_col())

write_xlsx(combine_table2, paste0(ABpath,"Output/table2.xlsx"))
