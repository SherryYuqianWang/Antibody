# Package ######################################################################
## Maths / Biostats -----------------------------------------------------------
library(deSolve) 
#library(MASS) #sample from multivariate normal distribution
library(finalfit) #table of regression result
library(survival) #coxph function
library(SurvMetrics)
## Plotting --------------------------------------------------------------------
library(ggplot2)
library(ggsci)
library(ggh4x)
#library(cowplot)
## Data processing -------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(writexl)
#library(purrr) #Functional Programming


rm(list=ls())
# Survival analysis ############################################################
## Data prep -------------------------------------------------------------------
#ba1_iga_cat <- read.csv("C:/Users/yuqian.wang/NTU_Sherry/4antibody/2data/ba1_iga_new.csv")
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

#full_ba1_iga <- full_ba1_iga %>% mutate(case_cat = cases/5920000*100000)

test <- full_ba1_iga %>% subset(doi_new==tstop & event==1) %>% arrange(cases)
summary(test$cases)

## Coxph model -----------------------------------------------------------------
ba1_iga_cont <- coxph(Surv(tstart, tstop, event) ~ antibody + log2(cases) , data = full_ba1_iga)

ba1_iga_cont <- coxph(Surv(tstart, tstop, event) ~ day28_ab + log2(cases) , data = full_ba1_iga)
                      #+ Gender + age_cat + Booster
#ba1_iga_cont <- coxph(Surv(tstart, tstop, event) ~ day28_ab + Gender + age_cat + Booster + log2(cases) , data = full_ba1_iga)
Cindex(coxph(Surv(tstart, tstop, event) ~ antibody + log2(cases) , data = full_ba1_iga), predicted)
             + Gender + age_cat + Booster 


bhest <- basehaz(ba1_iga_cont,centered = FALSE)
ba1_iga_cont$means #mean cases and mean antibody

#+ Gender + age_cat + Booster, data = full_ba1_iga ) 

summary(ba1_iga_cont) # higher ba1 IgA => lower risk of infection


mean_ab<-ba1_iga_cat %>% 
  group_by(time) %>% 
  summarise(mean_antibody = mean(antibody_new, na.rm = TRUE))


## Prediction ------------------------------------------------------------------
#rep_n <- length(seq(1,51,5))
rep_n=1
Tmin=1
#age_cat <- c("<60",">=60") # <60 and >=60
#booster <- c(0,1) # Moderna and Pfizer
#gender <- c("Female","Male") #female and male
#booster_int <- c(0,1,2) # 6months and 9months
cases <- c(2650,4350,8100) #low:600 cases(~10/100,000 people/day) #intermediate:5000 cases(~80/100,000 people/day) #high:12000 cases(~200/100,000 people/day)
tstop <- seq(9,339,10) #plot group by antibody
tstop <- seq(69,339,90) # plot group by time +21 days

pred_df<-data.frame()

#for (a in 1:2){
 #for (b in 1:2){
  for (c in 1:3){
   #for (g in 1:2){
     for (t in 1:length(tstop)){
       #for (i in 1:3){
  test.dat <- data.frame(#day28_ab = rep(seq(10,50,10),1),
                         antibody = rep(seq(1,51,2),1),#plot group by time
                         #antibody = rep(seq(10,50,10),1), #plot group by antibody
                         tstart = rep(1,rep_n),
                         tstop = rep(tstop[t],rep_n),
                         event = rep(0,rep_n),
                         #age_cat = rep(age_cat[a],rep_n),
                         #Booster = rep(booster[b],rep_n),
                         cases = rep(cases[c],rep_n),
                         #Gender = rep(gender[g],rep_n),
                         #booster_int = rep(booster_int[i],rep_n),
                         prob = rep(NA,rep_n))
  pred_df <- rbind(pred_df,test.dat)
    #   }
   #  }
   # }
 #  }
  }
}

#time-independent
preds <- predict(ba1_iga_cont, newdata = pred_df, type = "survival",se.fit = TRUE)

pred_df$prob <-preds$fit
pred_df$upr <- preds$fit + (1.96 * preds$se.fit)
pred_df$lwr <- preds$fit - (1.96 * preds$se.fit)


#time-varying
for (i in 1:nrow(pred_df)){
  
  #gender_ind <- ifelse(pred_df$Gender[i]=="Male",1,0)
  #age_ind <- ifelse(pred_df$age_cat[i]==">=60",1,0)
  #booster_ind <- ifelse(pred_df$Booster[i]=="Pfizer",1,0)
  
  #boost_interval <- pred_df$booster_int[i]
  tstop_ind <- pred_df$tstop[i]
  
  
    pars <- c(A0 <- pred_df$antibody[i],
              #k1 <- 0.13*exp(0.09),
              k2 <- 0.00614*(1+exp(0.2))/2)
              #k3 <- 0.07*exp(age_ind*0.21-gender_ind*0.1+pred_df$booster_int[i]*0.71))

  ind_ab <- ode_ab_surv(pars,tstop_ind)
  #ind_ab <- ode_ab_boost(pars,boost_interval,tstop_ind)
  #ind_ab_bind<-cbind(ind_ab1,ind_ab2[2],ind_ab3[2])
  #colnames(ind_ab_bind) <- c("time","3 months","6 months","9 months")
  #plot(plt)
  

  #plt <- ind_ab_bind %>% 
    #mutate(bashaz = c(0,diff(bhest[, 1])),
    #       cases = log2(ba1_iga_cat$cases[2:339])/1000000) %>%
  #  pivot_longer(!time, 
  #               names_to = 'Type', 
  #               values_to = 'Number')
  
  #ggplot(data=plt)+
  #  geom_line(aes(x = time,y = Number,colour=Type))+
    #xlab("Omicron IgA binding (%)")+
    #ylab("Protection against infection (%)")+
  #  theme_classic()
  
  #tstop_ind<-60
  
  pred_ind <- data.frame()
  pred_ind <- data.frame(#day28_ab = rep(seq(0,49,5),1),
    antibody = ind_ab[(1:tstop_ind),2],
    tstart = seq(Tmin,(tstop_ind),1),
    tstop = seq((Tmin+1),(tstop_ind+1),1),
    event=rep(0,(tstop_ind)),
    #age_cat=rep(pred_df$age_cat[i],(tstop_ind-21)),
    #Booster=rep(pred_df$Booster[i],(tstop_ind-21)),
    cases=rep(pred_df$cases[i],tstop_ind))
    #Gender=rep(pred_df$Gender[i],(tstop_ind-21)))
  
  preds <- predict(ba1_iga_cont, newdata = pred_ind, type = "survival",se.fit = TRUE)
  pred_df$prob[i] <- exp(sum(log(preds$fit)))
}

pred_plt <- pred_df %>% 
  mutate(#Booster = case_when(Booster == "0" ~ "Moderna",
         #                    Booster == "1" ~ "Pfizer"),
         #age_cat = case_when(age_cat == "<60" ~ "<60 yr",
         #                    age_cat == ">=60" ~ ">=60 yr"),
         #group = paste(Gender, age_cat, sep=", "),
         cases = case_when(cases == 2650 ~ "Low incidence",
                           cases == 4350 ~ "Medium incidence",
                           cases == 8100 ~ "High incidence"),
       #antibody = as.factor(antibody)) #plot group by antibody   
         #upr = ifelse(upr<=1,upr,1), #day28
         #lwr = ifelse(lwr<=0,0,lwr), #day28
         #day28_ab = as.factor(day28_ab)) #day28
        Duration = tstop+21,
        Duration = as.factor(Duration)) #plot group time
      
         #booster_int = case_when(booster_int == "0" ~ "6 months",
         #                        booster_int == "1" ~ "9 months",
         #                        booster_int == "2" ~ "Without fourth dose",),
         #tstop = case_when(tstop == "69" ~ 90,
         #                  tstop == "159" ~ 180,
         #                  tstop == "249" ~ 270,
         #                  tstop == "339" ~ 360)) %>%
#pivot_wider(names_from = booster_int, values_from = prob, values_fill = 0)
pred_plt$cases <- factor(pred_plt$cases, levels = c("Low incidence", "Medium incidence", "High incidence"))

#write_xlsx(pred_plt, "./output/prediction_result.xlsx")








for (i in 1:nrow(pred_df)){
  
  gender_ind <- 0
  age_ind <- 0
  #booster_ind <- ifelse(pred_df$Booster[i]=="Pfizer",1,0)
  
  boost_interval <- booster_int[i]
  tstop_ind <- 339
  
  
  pars <- c(A0 <- 1.07,
            k1 <- 0.13,
            k2 <- 0.00614*exp(0.2),
            k3 <- 0.07*exp(age_ind*0.21-gender_ind*0.1+booster_int[i]*0.71))
  
  #ind_ab <- ode_ab_surv(pars,tstop_ind)
  ind_ab <- ode_ab_boost(pars,boost_interval,tstop_ind)
  
  
  
  
  
  



bhest_plt <- bhest %>% 
  mutate(bashaz = c(0,diff(bhest[, 1])),
         cases = log2(ba1_iga_cat$cases[2:339])/1000000) %>%
  pivot_longer(cols = c('bashaz','cases'), 
               names_to = 'Type', 
               values_to = 'Number')

ggplot(data=bhest_plt)+
  geom_line(aes(x = time,y = Number,colour=Type))+
  #xlab("Omicron IgA binding (%)")+
  #ylab("Protection against infection (%)")+
  theme_classic()



pred_df$prob <-preds$fit
pred_df$upr <- preds$fit + (1.96 * preds$se.fit)
pred_df$lwr <- preds$fit - (1.96 * preds$se.fit)
#pred_df$cases <- as.factor(pred_df$cases)

pred_plt <- pred_df %>% 
  mutate(Booster = case_when(Booster == "0" ~ "Moderna",
                        Booster == "1" ~ "Pfizer"),
    age_cat = case_when(age_cat == "<60" ~ "<60 yr",
                       age_cat == ">=60" ~ ">=60 yr"),
    group = paste(Gender, age_cat, sep=", "),
    cases = case_when(cases == 600 ~ "Low transmission",
                      cases == 5000 ~ "Intermediate transmission",
                      cases == 12000 ~ "High transmission"))
    #upr = ifelse(upr<=1,upr,1),
    #lwr = ifelse(lwr<=0,0,lwr))



## Figure ----------------------------------------------------------------------
#antibody level(time varying)
low_20 <- 11 
mod_20 <- 9
high_20 <- 23

low_40 <- 35
mod_40 <- 25
high_40 <- 43

#antibody peak (time varying)
low_20 <- 215 
mod_20 <- 155
high_20 <- 95

low_40 <- 290
mod_40 <- 225
high_40 <- 170


#antibody peak (day28
low_20 <- 225 #+10
mod_20 <- 165  #+10
high_20 <- 105  #+5

low_40 <- 365  #+70
mod_40 <- 280  #+55
high_40 <- 190 #+20

#antibody=40
#line_dat = data.frame(cases=c("Low incidence","Medium incidence","High incidence"),
#                      xv=c(low_40,mod_40, high_40),
#                      xendv=c(low_40,mod_40,high_40),
#                      yv=c(0,0, 0),
#                      yendv=c(80,80,80),
#                      xh=c(0,0,0),
#                      xendh=c(low_40,mod_40, high_40),
#                      yh=c(80,80, 80),
#                      yendh=c(80,80,80))

#line_dat$cases <- factor(line_dat$cases, levels = c("Low incidence", "Medium incidence", "High incidence"))

#antibody=20 and antibody=40
line_dat = data.frame(cases=c("Low incidence","Medium incidence","High incidence","Low incidence","Medium incidence","High incidence"),
                      xv=c(low_40,mod_40, high_40,low_20,mod_20, high_20),
                      xendv=c(low_40,mod_40,high_40,low_20,mod_20, high_20),
                      yv=c(0,0, 0),
                      yendv=c(80,80,80),
                      xh=c(0,0,0),
                      xendh=c(low_40,mod_40, high_40,low_20,mod_20, high_20),
                      yh=c(80,80, 80),
                      yendh=c(80,80,80))

line_dat$cases <- factor(line_dat$cases, levels = c("Low incidence", "Medium incidence", "High incidence"))






ggplot(data=pred_plt)+
  geom_smooth(aes(x = antibody, y = prob*100,linetype=Duration),linewidth=0.7)+
  #geom_line(aes(x = antibody, y = prob*100,colour=group,linetype=Booster))+
  #geom_ribbon(aes(x = day28_ab, ymin = lwr*100, ymax = upr*100),alpha=0.5,fill="#003399FF") +
  #labs(colour="Group") +
  geom_segment(data=line_dat, aes(x = xh, xend = xendh, y = yh, yend = yendh), linetype="dashed", color = "darkblue", linewidth=0.3)+
  geom_segment(data=line_dat, aes(x = xv, xend = xendv, y = yv, yend = yendv), linetype="dashed", color = "darkblue", linewidth=0.3)+
  #scale_color_manual(values=c("#0099CCFF","#003399FF","#EE4C97FF","#A20056FF"))+
  facet_grid2(~cases,axes = "all")+
  #xlim(0,51)+
  scale_linetype_manual(values=c("solid","dotted","dotdash","longdash")) +  
  scale_y_continuous(breaks=seq(0,100,by=10),labels = expression(0,10,20,30,40,50,60,70,80,90,100),limits=c(0,100)) +
  xlab("BA.1 IgA binding (%)")+
  ylab("Protection against infection (%)")+
  labs(linetype = "Time after booster vaccination (days)")+
  theme_classic()+
  theme(legend.position="bottom")
#"dashed",
ggsave("plot/protection_updated_ab.png", width = 12, height = 4)


ggplot(data=pred_plt)+
  geom_smooth(aes(x = (tstop+21), y = prob*100 , linetype = antibody),linewidth=0.7)+
  #geom_smooth(aes(x = (tstop+28), y = prob*100 , linetype = day28_ab),linewidth=0.7)+
  geom_segment(data=line_dat, aes(x = xh, xend = xendh, y = yh, yend = yendh), linetype="dashed", color = "darkblue", linewidth=0.3)+
  geom_segment(data=line_dat, aes(x = xv, xend = xendv, y = yv, yend = yendv), linetype="dashed", color = "darkblue", linewidth=0.3)+
  #geom_text(data = line_dat, aes(x = xv, y = yv, label = xv), vjust = 0, hjust = 0) +
  #scale_x_continuous(breaks = c(100, 200, 300,xv)) +
  #geom_vline(xintercept = c(50,60))+
  scale_y_continuous(breaks=seq(0,100,by=10),labels = expression(0,10,20,30,40,50,60,70,80,90,100),limits=c(0.0,100)) +
  scale_linetype_manual(values=c("solid","dotted","dashed","dotdash","longdash")) +  
  facet_grid2(~cases,axes = "all",scales = "free_x")+
  facetted_pos_scales(
    x = list(
      cases == "Low incidence" ~ scale_x_continuous(breaks = c(0,50,100,150,350,low_40,low_20)),
      cases == "Medium incidence" ~ scale_x_continuous(breaks = c(0,50,100,300,350,mod_40,mod_20)),
      #cases == "Low incidence" ~ scale_x_continuous(breaks = c(0,50,100,150,low_40,low_20)), #day28
      #cases == "Medium incidence" ~ scale_x_continuous(breaks = c(0,50,100,350,mod_40,mod_20)), #day28
      cases == "High incidence" ~ scale_x_continuous(breaks = c(0,50,250,300,350,high_40,high_20))))+
  xlab("Time after booster vaccination (days)")+
  ylab("Protection against infection (%)")+
  labs(linetype = "BA.1 IgA binding (%)")+
  theme_classic()+
  theme(legend.position="bottom")


ggsave("plot/protection_antibody.png", width = 12, height = 4)

