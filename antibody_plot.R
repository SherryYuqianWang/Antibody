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
setwd("C:/Users/yuqian.wang/NTU_Sherry/4antibody/6code")
source("antibody_function.R") #function for antibody plots

# Setting ######################################################################
Tmin <- 0
Tmax <- 360

step_size <- 1
times <- c(seq(Tmin,Tmax,step_size))

# Population fit ###############################################################
## Import data -----------------------------------------------------------------
pop_fit <- list()
pop_fit[[1]] <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/4antibody/3monolix/r11_wt_igg_03/populationParameters.txt", row.names = 1)
pop_fit[[2]] <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/4antibody/3monolix/r11_ba1_igg_14/populationParameters.txt", row.names = 1)
pop_fit[[3]] <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/4antibody/3monolix/r11_wt_iga_07/populationParameters.txt", row.names = 1)
pop_fit[[4]] <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/4antibody/3monolix/r11_ba1_iga_03/populationParameters.txt", row.names = 1)

## Generate figure -------------------------------------------------------------
num = 500
ab_list <- c("IgG binding (%)","IgG binding (%)","IgA binding (%)","IgA binding (%)")
group_df <- list()

# age 
group_df[[1]] <- data.frame(group = c("Pf booster, early infected, age <60, male",
                                      "Pf booster, late infected, age <60, male",
                                      "Pf booster, non infected, age <60, male",
                                      "Pf booster, early infected, age >=60, male",
                                      "Pf booster, late infected, age >=60, male",
                                      "Pf booster, non infected, age >=60, male"),
                            booster = rep(1,6),
                            infect_late = c(0, 1, 0, 0, 1, 0),
                            infect_un = c(0, 0, 1, 0, 0, 1),
                            infect = c(0, 1, 2, 0, 1, 2),
                            gender = rep(0,6),
                           age = c(0,0,0,1,1,1))
# gender 
group_df[[2]] <- data.frame(group = c("Pf booster, early infected, age <60, male",
                                      "Pf booster, late infected, age <60, male",
                                      "Pf booster, non infected, age <60, male",
                                      "Pf booster, early infected, age <60, female",
                                      "Pf booster, late infected, age <60, female",
                                      "Pf booster, non infected, age <60, female"),
                            booster = rep(1,6),
                            infect_late = c(0, 1, 0, 0, 1, 0),
                            infect_un = c(0, 0, 1, 0, 0, 1),
                            infect = c(0, 1, 2, 0, 1, 2),
                            gender = c(1, 1, 1, 0, 0, 0),
                            age = rep(0,6))
# booster
group_df[[3]] <- data.frame(group = c("Pf booster, early infected, male",
                                     "Pf booster, late infected, male",
                                     "Pf booster, non infected, male",
                                     "Mo booster, early infected, male",
                                     "Mo booster, late infected, male",
                                     "Mo booster, non infected, male"),
                           booster = c(1, 1, 1, 0, 0, 0),
                           infect_late = c(0, 1, 0, 0, 1, 0),
                           infect_un = c(0, 0, 1, 0, 0, 1),
                           infect = c(0, 1, 2, 0, 1, 2),
                           gender = rep(0,6),
                           age = rep(0,6))


plot_list = list()

for (i in 1:4){
  pop <- pop_fit[[i]]
  ab <- ab_list[i]

    for (j in 1:3){
    group_list <- split(group_df[[j]],seq(nrow(group_df[[j]]))) 
    
    total_AB <- map(group_list,simulate_AB)
    
    test_AB <- list_rbind(total_AB) %>% mutate(Booster =
                                                 case_when(Booster == 0 ~ "Moderna", 
                                                           Booster == 1 ~ "Pfizer"), 
                                               Infect = 
                                                 case_when(Infect == 0 ~ "Early infection",
                                                           Infect == 1 ~ "Late infection",
                                                           Infect == 2 ~ "Uninfected"),
                                               Age = 
                                                 case_when(Age == 0 ~ "<60 yr",
                                                           Age == 1 ~ ">=60 yr"),
                                               Gender = 
                                                 case_when(Gender == 0 ~ "Female", 
                                                           Gender == 1 ~ "Male"))
    if (j==1 & ab=="IgG binding (%)"){
      p <- ggplot() +
      geom_line(data=test_AB,aes(x=time,y=best_ab,color = Age),lwd=1) +
      geom_ribbon(data=test_AB,aes(x=time,ymin=q1_ab,ymax=q3_ab, fill = Age),alpha=0.15) +
      geom_ribbon(data=test_AB,aes(x=time,ymin=low_ab,ymax=high_ab, fill = Age),alpha=0.05) +
      xlab("Time after booster") +
      ylab(ab_list[i])  +
      scale_y_continuous(breaks=seq(10,70,by=10),labels = expression(10,20,30,40,50,60,70),limits=c(0,70)) +
      scale_color_manual(values=c("#a2e7f3","#00429d"))+
      scale_fill_manual(values=c("#a2e7f3","#00429d"))+
      facet_grid(. ~ Infect)+
      theme_classic()+
      theme(axis.title.x=element_blank(),
            legend.position="none")
      
    } else if (j==1 & ab!="IgG binding (%)"){
      p <- ggplot() +
        geom_line(data=test_AB,aes(x=time,y=best_ab,color = Age),lwd=1) +
        geom_ribbon(data=test_AB,aes(x=time,ymin=q1_ab,ymax=q3_ab, fill = Age),alpha=0.15) +
        geom_ribbon(data=test_AB,aes(x=time,ymin=low_ab,ymax=high_ab, fill = Age),alpha=0.05) +
        xlab("Time after booster") +
        ylab(ab_list[i])  +
        scale_y_continuous(breaks=seq(10,70,by=10),labels = expression(10,20,30,40,50,60,70),limits=c(0,70)) +
        scale_color_manual(values=c("#a2e7f3","#00429d"))+
        scale_fill_manual(values=c("#a2e7f3","#00429d"))+
        facet_grid(. ~ Infect)+
        theme_classic()+
        theme(axis.title.x=element_blank())
      
    } else if (j==2 & ab=="IgG binding (%)"){
      p <- ggplot() +
      geom_line(data=test_AB,aes(x=time,y=best_ab,color = Gender),lwd=1) +
      geom_ribbon(data=test_AB,aes(x=time,ymin=q1_ab,ymax=q3_ab, fill = Gender),alpha=0.2) +
      geom_ribbon(data=test_AB,aes(x=time,ymin=low_ab,ymax=high_ab, fill = Gender),alpha=0.08) +
      xlab("Time after booster") +
      ylab(ab_list[i])  +
      scale_y_continuous(breaks=seq(10,70,by=10),labels = expression(10,20,30,40,50,60,70),limits=c(0,70)) +
      scale_color_lancet()+
      scale_fill_lancet()+
      facet_grid(. ~ Infect)+
      theme_classic()+
        theme(axis.title.x=element_blank(),
              legend.position="none")
      
    } else if (j==2 & ab!="IgG binding (%)"){
      p <- ggplot() +
        geom_line(data=test_AB,aes(x=time,y=best_ab,color = Gender),lwd=1) +
        geom_ribbon(data=test_AB,aes(x=time,ymin=q1_ab,ymax=q3_ab, fill = Gender),alpha=0.2) +
        geom_ribbon(data=test_AB,aes(x=time,ymin=low_ab,ymax=high_ab, fill = Gender),alpha=0.08) +
        xlab("Time after booster") +
        ylab(ab_list[i])  +
        scale_y_continuous(breaks=seq(10,70,by=10),labels = expression(10,20,30,40,50,60,70),limits=c(0,70)) +
        scale_color_lancet()+
        scale_fill_lancet()+
        facet_grid(. ~ Infect)+
        theme_classic()+
        theme(axis.title.x=element_blank())
      
      } else if (j==3 & (ab=="IgG binding (%)")){
      p <- ggplot() +
      geom_line(data=test_AB,aes(x=time,y=best_ab,color = Booster),lwd=1) +
      geom_ribbon(data=test_AB,aes(x=time,ymin=q1_ab,ymax=q3_ab, fill = Booster),alpha=0.2) +
      geom_ribbon(data=test_AB,aes(x=time,ymin=low_ab,ymax=high_ab, fill = Booster),alpha=0.08) +
      xlab("Time after booster") +
      ylab(ab_list[i])  +
      scale_y_continuous(breaks=seq(10,70,by=10),labels = expression(10,20,30,40,50,60,70),limits=c(0,70)) +
      facet_grid(. ~ Infect)+
      theme_classic()+
        theme(legend.position="none")
  
      } else {
        p <- ggplot() +
          geom_line(data=test_AB,aes(x=time,y=best_ab,color = Booster),lwd=1) +
          geom_ribbon(data=test_AB,aes(x=time,ymin=q1_ab,ymax=q3_ab, fill = Booster),alpha=0.2) +
          geom_ribbon(data=test_AB,aes(x=time,ymin=low_ab,ymax=high_ab, fill = Booster),alpha=0.08) +
          xlab("Time after booster") +
          ylab(ab_list[i])  +
          scale_y_continuous(breaks=seq(10,70,by=10),labels = expression(10,20,30,40,50,60,70),limits=c(0,70)) +
          facet_grid(. ~ Infect)+
          theme_classic()
    }
    plot_list[[(i-1)*3+j]] = p
    
  }
}

  
ggdraw() +
  draw_plot(plot_list[[3]], x = 0.02, y = 0.72, width = 0.4, height = 0.24) +
  draw_plot(plot_list[[4]], x = 0.02, y = 0.45, width = 0.4, height = 0.22) +
  draw_plot(plot_list[[5]], x = 0.02, y = 0.24, width = 0.4, height = 0.22) +
  draw_plot(plot_list[[6]], x = 0.02, y = 0, width = 0.4, height = 0.24) +
  draw_plot(plot_list[[9]], x = 0.47, y = 0.72, width = 0.53, height = 0.24) +
  draw_plot(plot_list[[10]], x = 0.47, y = 0.45, width = 0.53, height = 0.22) +
  draw_plot(plot_list[[11]], x = 0.47, y = 0.24, width = 0.53, height = 0.22) +
  draw_plot(plot_list[[12]], x = 0.47, y = 0, width = 0.53, height = 0.24) +
  draw_plot_label(label = c("WT IgG", "Omicron IgG","WT IgA", "Omicron IgA" ), size = 14,
                  x = c(0.037,0.0,0.487,0.45), y = c(rep(c(0.995, 0.705),2)))

ggsave("plot/population2.png", width = 9.5, height = 7.5,bg = "white")

# Individual fit ###############################################################
## Import data -----------------------------------------------------------------
original <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/4antibody/2data/P360_Sflow_For Keisuke_17Jan2024_remove74.csv")
original$booster_trans <- strptime(as.character(original$Booster.Vaccination.date), "%d/%m/%Y")
original_cov <- list()
original_cov[[1]] <- original %>% select(Code, booster_trans, DoI, WT.IgG) %>% group_by(Code) %>% slice(2)
original_cov[[2]] <- original %>% select(Code, booster_trans, DoI, BA1.IgG) %>% group_by(Code) %>% slice(2)
original_cov[[3]] <- original %>% select(Code, booster_trans, DoI, WT.IgA) %>% group_by(Code) %>% slice(2)
original_cov[[4]] <- original %>% select(Code, booster_trans, DoI, BA1.IgA) %>% group_by(Code) %>% slice(2)

sg_covid <- list()
for (i in 1:4){
  sg_covid[[i]] <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/4antibody/2data/weekly_cases_per_million.csv")
  sg_covid[[i]]$date_trans <- strptime(as.character(sg_covid[[i]]$date), "%Y-%m-%d")
}

Estimated <- list()
Estimated[[1]] <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/4antibody/3monolix/r11_wt_igg_03/IndividualParameters/estimatedIndividualParameters.txt", sep = ",", comment.char = "", header = T)
Estimated[[2]] <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/4antibody/3monolix/r11_ba1_igg_14/IndividualParameters/estimatedIndividualParameters.txt", sep = ",", comment.char = "", header = T)
Estimated[[3]] <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/4antibody/3monolix/r11_wt_iga_07/IndividualParameters/estimatedIndividualParameters.txt", sep = ",", comment.char = "", header = T)
Estimated[[4]] <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/4antibody/3monolix/r11_ba1_iga_03/IndividualParameters/estimatedIndividualParameters.txt", sep = ",", comment.char = "", header = T)

Simulated <- list()
Simulated[[1]] <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/4antibody/3monolix/r11_wt_igg_03/IndividualParameters/simulatedIndividualParameters.txt", sep = ",", comment.char = "", header = T)
Simulated[[2]] <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/4antibody/3monolix/r11_ba1_igg_14/IndividualParameters/simulatedIndividualParameters.txt", sep = ",", comment.char = "", header = T)
Simulated[[3]] <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/4antibody/3monolix/r11_wt_iga_07/IndividualParameters/simulatedIndividualParameters.txt", sep = ",", comment.char = "", header = T)
Simulated[[4]] <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/4antibody/3monolix/r11_ba1_iga_03/IndividualParameters/simulatedIndividualParameters.txt", sep = ",", comment.char = "", header = T)

## Tables ----------------------------------------------------------------------

#Generate dataset with k1,k2,k3,k4
ind_AB <- list()
ind_AB <- pmap(list(Estimated,original_cov,sg_covid),ind_ab_long)

#Generate dataset with k1,k2
ind_AB <- list()
ind_AB <- map2(Estimated,original_cov,ind_ab_long_cluster)


sheets <- list("wt_igg" = ind_AB[[1]], "ba1_igg" = ind_AB[[2]], "wt_iga" = ind_AB[[3]], "ba1_iga" = ind_AB[[4]])
write_xlsx(sheets, "./output/antibody_individual_protection_80%.xlsx")

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

write_xlsx(combine_reg, "./output/regression_result.xlsx")


# Survival analysis ############################################################
## Data prep -------------------------------------------------------------------
ba1_iga_cat <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/4antibody/2data/ba1_iga_new.csv")

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
rep_n <- length(seq(0,49,5))

#age_cat <- c(">=60","<60")
#booster <- c("1","0")
#gender <- c("Male","Female")
cases <- c(600,5000,12000) #low:600 cases(~10/100,000 people/day) #intermediate:5000 cases(~80/100,000 people/day) #high:12000 cases(~200/100,000 people/day)



pred_df<-data.frame()

#for (a in 1:2){
#for (t in 1:2){
for (c in 1:3){
  # for (g in 1:2){
  test.dat <- data.frame(antibody = rep(seq(0,49,5),1),
                         tstart=rep(1,rep_n),
                         tstop=rep(91,rep_n),
                         event=rep(0,rep_n),
                         # age_cat=rep(age_cat[a],rep_n),
                         # Booster=rep(booster[b],rep_n),
                         cases=rep(cases[c],rep_n))
  # Gender=rep(gender[g],rep_n))
  pred_df <- rbind(pred_df,test.dat)
  #}
  #}
  #}
}


preds <- predict(ba1_iga_cont, newdata = pred_df, type = "survival",se.fit = TRUE)

pred_df$prob <-preds$fit
pred_df$upr <- preds$fit + (1.96 * preds$se.fit)
pred_df$lwr <- preds$fit - (1.96 * preds$se.fit)
#pred_df$cases <- as.factor(pred_df$cases)

pred_df <- pred_df %>% 
  mutate(#Booster = case_when(Booster == "0" ~ "Moderna",
    #                    Booster == "1" ~ "Pfizer"),
    #age_cat = case_when(age_cat == "<60" ~ "<60 yr",
    #                   age_cat == ">=60" ~ ">=60 yr"),
    cases = case_when(cases == 600 ~ "Low transmission",
                      cases == 5000 ~ "Intermediate transmission",
                      cases == 12000 ~ "High transmission"),
    upr = ifelse(upr<=1,upr,1),
    lwr = ifelse(lwr<=0,0,lwr))

pred_df$cases <- factor(pred_df$cases, levels = c("Low transmission", "Intermediate transmission", "High transmission"))


## Figure ----------------------------------------------------------------------
#ba1_IgA_3 months
mod_80 <- 10
high_80 <- 23.6

#ba1_IgA_6 months
mod_80 <- 20
high_80 <- 34

line_dat = data.frame(cases=c("Low transmission","Intermediate transmission","High transmission"),
                      xv=c(-1,mod_80, high_80),
                      xendv=c(-1,mod_80,high_80),
                      yv=c(-1,0, 0),
                      yendv=c(-1,80, 80),
                      xh=c(-1,0,0),
                      xendh=c(-1,mod_80, high_80),
                      yh=c(-1,80, 80),
                      yendh=c(-1,80,80))

line_dat$cases <- factor(line_dat$cases, levels = c("Low transmission", "Intermediate transmission", "High transmission"))


ggplot(data=pred_df)+
  geom_line(aes(x = antibody,y = prob*100))+
  geom_ribbon(aes(x = antibody, ymin = lwr*100, ymax = upr*100),alpha=0.5,fill="#003399FF") +
  labs(colour="Group") +
  geom_segment(data=line_dat, aes(x = xh, xend = xendh, y = yh, yend = yendh), linetype="dashed", color = "darkblue", linewidth=0.3)+
  geom_segment(data=line_dat, aes(x = xv, xend = xendv, y = yv, yend = yendv), linetype="dashed", color = "darkblue", linewidth=0.3)+
  facet_grid2(~cases,axes = "all")+
  ylim(0,100)+
  xlim(0,49)+
  scale_y_continuous(breaks=seq(0,100,by=10),labels = expression(0,10,20,30,40,50,60,70,80,90,100),limits=c(0.0,100)) +
  xlab("Omicron IgA binding (%)")+
  ylab("Protection against infection (%)")+
  theme_classic()

ggsave("plot/protect_ba1_3months.png", width = 8, height = 3)




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

write_xlsx(combine_table2, "./output/table2.xlsx")
