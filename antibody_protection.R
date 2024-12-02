# load packages ----------------------------------------------------------------



# Data prep ---------------------------------------------------------------
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

# Survival for continuous antibody levels -------------------------------
ba1_iga_cont <- coxph(Surv(tstart, tstop, event) ~ antibody + 
                        log2(cases) , data = full_ba1_iga)
                       #+ Gender + age_cat + Booster, data = full_ba1_iga ) 

summary(ba1_iga_cont) # higher ba1 IgA => lower risk of infection


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
   # }
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
