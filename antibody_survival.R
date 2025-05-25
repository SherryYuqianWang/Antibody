# load packages ----------------------------------------------------------------
library(broom)             # helps to tidy up model outputs
library(survival)          # coxph function
library(coxme)             # cox ph with mixed effects
library(tidyverse)
library(kableExtra)
library(timereg)           # to test for proportionality
library(survminer)         # for kaplan meier
library(rstpm2)            # tvc function
library(rms)               # restricted cubic splines
library(randomForestSRC)   # survival through random survival forests
library(ggplot2)
library(forestplot)
library(coxphf)            # survival model with Firth's penalized likelihood
library(survminer)         # for cut points
library(BestSurvivalCuts)  # cut points based on AIC
library(kml)               # longitudinal k-means
library(clValid)           # test cluster score
library(corrr)             # correlation matrix
library(ggcorrplot)
library(mgcv)              # for GAM
library(pROC)
library(autothresholdr)    # data-driven threshold estimates
library(gridExtra)         # to arrange ggplots
library(survcompare)       # 
library(timeROC)           # time varying ROC

# Data prep ---------------------------------------------------------------
# Read in original data files
wt_igg_cat <- read.csv("~/Dropbox (Personal)/Ejima-lab (LKC)/Projects/4. Antibody Modelling/wt_igg_new.csv")
wt_iga_cat <- read.csv("~/Dropbox (Personal)/Ejima-lab (LKC)/Projects/4. Antibody Modelling/wt_iga_new.csv")
ba1_igg_cat <- read.csv("~/Dropbox (Personal)/Ejima-lab (LKC)/Projects/4. Antibody Modelling/ba1_igg_new.csv")
ba1_iga_cat <- read.csv("~/Dropbox (Personal)/Ejima-lab (LKC)/Projects/4. Antibody Modelling/ba1_iga_new.csv")

# prep data frames for survival models
wt_igg_cat$doi_new[is.na(wt_igg_cat$doi_new)] <- 339
wt_igg_df_cat <- filter(wt_igg_cat, time <= doi_new)
wt_igg_df_cat$Infect[wt_igg_df_cat$Infect == 0] <- 1
wt_igg_df_cat$Infect[wt_igg_df_cat$Infect == 2] <- 0
#wt_igg_df_cat <- rename(wt_igg_df_cat, WT.IgG = antibody)

wt_iga_cat$doi_new[is.na(wt_iga_cat$doi_new)] <- 339
wt_iga_df_cat <- filter(wt_iga_cat, time <= doi_new)
wt_iga_df_cat$Infect[wt_iga_df_cat$Infect == 0] <- 1
wt_iga_df_cat$Infect[wt_iga_df_cat$Infect == 2] <- 0
#wt_iga_df_cat <- rename(wt_iga_df_cat, WT.IgA = antibody)

ba1_igg_cat$doi_new[is.na(ba1_igg_cat$doi_new)] <- 339
ba1_igg_df_cat <- filter(ba1_igg_cat, time <= doi_new)
ba1_igg_df_cat$Infect[ba1_igg_df_cat$Infect == 0] <- 1
ba1_igg_df_cat$Infect[ba1_igg_df_cat$Infect == 2] <- 0
#ba1_igg_df_cat <- rename(ba1_igg_df_cat, BA1.IgG = antibody)

ba1_iga_cat$doi_new[is.na(ba1_iga_cat$doi_new)] <- 339
ba1_iga_df_cat <- filter(ba1_iga_cat, time <= doi_new)
ba1_iga_df_cat$Infect[ba1_iga_df_cat$Infect == 0] <- 1
ba1_iga_df_cat$Infect[ba1_iga_df_cat$Infect == 2] <- 0
#ba1_iga_df_cat <- rename(ba1_iga_df_cat, BA1.IgA = antibody)

# separate time-dependent and time-independent variables
time_ind_wt_igg <- wt_igg_df_cat %>% dplyr::select(id, Gender, Infect, age_cat, Booster, 
                                            doi_new) %>% distinct()
time_ind_wt_igg <- tmerge(data1 = time_ind_wt_igg, data2 = time_ind_wt_igg, id = id, event = event(doi_new, Infect))
time_dep_wt_igg <- wt_igg_df_cat %>% dplyr::select(id, antibody, time, cases, 
                                            antibody_cat4, day28_cat4)
full_wt_igg <- tmerge(data1 = time_ind_wt_igg, data2 = time_dep_wt_igg, id = id, 
                      antibody = tdc(time, antibody), cases = tdc(time, cases), 
                      cat4 = tdc(time, antibody_cat4), day28_cat4 = tdc(time, day28_cat4))
full_wt_igg <- filter(full_wt_igg, tstart != 0)

time_ind_wt_iga <- wt_iga_df_cat %>% dplyr::select(id, Gender, Infect, age_cat, Booster, 
                                            doi_new) %>% distinct()
time_ind_wt_iga <- tmerge(data1 = time_ind_wt_iga, data2 = time_ind_wt_iga, id = id, event = event(doi_new, Infect))
time_dep_wt_iga <- wt_iga_df_cat %>% dplyr::select(id, antibody, time, cases, 
                                            antibody_cat4, day28_cat4)
full_wt_iga <- tmerge(data1 = time_ind_wt_iga, data2 = time_dep_wt_iga, id = id, 
                      antibody = tdc(time, antibody), cases = tdc(time, cases), 
                      cat4 = tdc(time, antibody_cat4), day28_cat4 = tdc(time, day28_cat4))
full_wt_iga <- filter(full_wt_iga, tstart != 0)

time_ind_ba1_igg <- ba1_igg_df_cat %>% dplyr::select(id, Gender, Infect, age_cat, Booster, 
                                              doi_new) %>% distinct()
time_ind_ba1_igg <- tmerge(data1 = time_ind_ba1_igg, data2 = time_ind_ba1_igg, id = id, event = event(doi_new, Infect))
time_dep_ba1_igg <- ba1_igg_df_cat %>% dplyr::select(id, antibody, time, cases, 
                                              antibody_cat4, day28_cat4)
full_ba1_igg <- tmerge(data1 = time_ind_ba1_igg, data2 = time_dep_ba1_igg, id = id, 
                       antibody = tdc(time, antibody), cases = tdc(time, cases), 
                       cat4 = tdc(time, antibody_cat4), day28_cat4 = tdc(time, day28_cat4))
full_ba1_igg <- filter(full_ba1_igg, tstart != 0)

time_ind_ba1_iga <- ba1_iga_df_cat %>% dplyr::select(id, Gender, Infect, age_cat, Booster, 
                                              doi_new) %>% distinct()
time_ind_ba1_iga <- tmerge(data1 = time_ind_ba1_iga, data2 = time_ind_ba1_iga, id = id, event = event(doi_new, Infect))
time_dep_ba1_iga <- ba1_iga_df_cat %>% dplyr::select(id, antibody, time, cases, 
                                              antibody_cat4, day28_cat4)
full_ba1_iga <- tmerge(data1 = time_ind_ba1_iga, data2 = time_dep_ba1_iga, id = id, 
                       antibody = tdc(time, antibody), cases = tdc(time, cases), 
                       cat4 = tdc(time, antibody_cat4), day28_cat4 = tdc(time, day28_cat4))
full_ba1_iga <- filter(full_ba1_iga, tstart != 0)

# day 28 categorical antibody data
d28_wt_igg_df <- wt_igg_df_cat %>% group_by(id) %>% slice(1) %>%
  mutate(day28_cat4 = recode(day28_cat4, "Low25%" = "Low", "High25%" = "High",
                             "25-75%" = "Medium")) %>%
  mutate(antibody_cat4 = recode(antibody_cat4, "Low25%" = "Low", "High25%" = "High",
                             "25-75%" = "Medium")) %>%
  rename(D28_antibody_level = day28_cat4, Antibody_cat = antibody_cat4)
d28_wt_igg_df$D28_antibody_level <- 
  factor(d28_wt_igg_df$D28_antibody_level, c("Medium", "High", "Low"))

d28_wt_iga_df <- wt_iga_df_cat %>% group_by(id) %>% slice(1) %>%
  mutate(day28_cat4 = recode(day28_cat4, "Low25%" = "Low", "High25%" = "High",
                             "25-75%" = "Medium")) %>%
  mutate(antibody_cat4 = recode(antibody_cat4, "Low25%" = "Low", "High25%" = "High",
                                "25-75%" = "Medium")) %>%
  rename(D28_antibody_level = day28_cat4, Antibody_cat = antibody_cat4)
d28_wt_iga_df$D28_antibody_level <- 
  factor(d28_wt_iga_df$D28_antibody_level, c("Medium", "High", "Low"))

d28_ba1_igg_df <- ba1_igg_df_cat %>% group_by(id) %>% slice(1) %>%
  mutate(day28_cat4 = recode(day28_cat4, "Low25%" = "Low", "High25%" = "High",
                             "25-75%" = "Medium")) %>%
  mutate(antibody_cat4 = recode(antibody_cat4, "Low25%" = "Low", "High25%" = "High",
                                "25-75%" = "Medium")) %>%
  rename(D28_antibody_level = day28_cat4, Antibody_cat = antibody_cat4)
d28_ba1_igg_df$D28_antibody_level <- 
  factor(d28_ba1_igg_df$D28_antibody_level, c("Medium", "High", "Low"))

d28_ba1_iga_df <- ba1_iga_df_cat %>% group_by(id) %>% slice(1) %>%
  mutate(day28_cat4 = recode(day28_cat4, "Low25%" = "Low", "High25%" = "High",
                             "25-75%" = "Medium")) %>%
  mutate(antibody_cat4 = recode(antibody_cat4, "Low25%" = "Low", "High25%" = "High",
                                "25-75%" = "Medium")) %>%
  rename(D28_antibody_level = day28_cat4, Antibody_cat = antibody_cat4)
d28_ba1_iga_df$D28_antibody_level <- 
  factor(d28_ba1_iga_df$D28_antibody_level, c("Low", "Medium", "High"))

# day 28 raw data
ab_df <- read.csv("~/Dropbox (Personal)/Ejima-lab (LKC)/Projects/4. Antibody Modelling/Data/Original/P360_Sflow_For Keisuke_17Jan2024_remove74.csv")
d28_ab_df <- filter(ab_df, Day == 28) %>% cbind(cases = d28_ba1_igg_df$cases)
d28_ab_df$DoI[is.na(d28_ab_df$DoI)] <- 365

# day 28 with time-updated cases
time_ind_d28_ab <- d28_ab_df %>% dplyr::select(Code, Gender, X1st.Booster.vaccine.type, 
                                        WT.IgG, WT.IgA, BA1.IgG, BA1.IgA, DoI,
                                        inf, age_cat)
time_ind_d28_ab <- tmerge(data1 = time_ind_d28_ab, data2 = time_ind_d28_ab, id = Code, event = event(DoI, inf))
time_dep_d28_ab <- wt_iga_df_cat %>% dplyr::select(id, time, cases)
full_d28_ab <- tmerge(data1 = time_ind_d28_ab, data2 = time_dep_d28_ab, id = id, 
                      cases = tdc(time, cases))
full_d28_ab <- filter(full_d28_ab, tstart != 0)

# Day 28 survival (continuous) -------------------------------------
# antibody levels
wt_igg_d28 <- coxph(Surv(tstart, tstop, event) ~ WT.IgG + Gender + age_cat 
                    + X1st.Booster.vaccine.type + log2(cases), data = full_d28_ab)
wt_iga_d28 <- coxph(Surv(tstart, tstop, event) ~ WT.IgA + Gender + age_cat 
                    + X1st.Booster.vaccine.type + log2(cases), data = full_d28_ab)
ba1_igg_d28 <- coxph(Surv(tstart, tstop, event) ~ BA1.IgG + Gender + age_cat 
                     + X1st.Booster.vaccine.type + log2(cases), data = full_d28_ab)
ba1_iga_d28 <- coxph(Surv(tstart, tstop, event) ~ BA1.IgA + Gender + age_cat 
                     + X1st.Booster.vaccine.type + log2(cases), data = full_d28_ab)
summary(wt_igg_d28)
summary(wt_iga_d28)
summary(ba1_igg_d28)
summary(ba1_iga_d28)

# Other covariate survival (univariate and multivariate adjusting for incidence) ------------------------
all_hr <- coxph(Surv(tstart, tstop, event) ~ WT.IgG + WT.IgA + BA1.IgG + BA1.IgA + 
                  Gender + age_cat + X1st.Booster.vaccine.type + log2(cases), 
                data = full_d28_ab)
summary(all_hr)

gender  <- coxph(Surv(tstart, tstop, event) ~ Gender, data = full_d28_ab)
age     <- coxph(Surv(tstart, tstop, event) ~ age_cat, data = full_d28_ab)
booster <- coxph(Surv(tstart, tstop, event) ~ X1st.Booster.vaccine.type, data = full_d28_ab)
summary(gender)
summary(age)
summary(booster)

# Tests for proportionality --------------------------------------------
zph_wt_igg <- cox.zph(wt_igg_d28, transform = "identity")
zph_wt_igg
plot(zph_wt_igg[1])
abline(h=wt_igg_d28$coefficients[1], col=3, lwd=2, lty=2)

zph_wt_iga <- cox.zph(wt_iga_d28)
zph_wt_iga
plot(zph_wt_iga[1])
abline(h=wt_iga_d28$coefficients[1], col=3, lwd=2, lty=2)

zph_ba1_igg <- cox.zph(ba1_igg_d28)
zph_ba1_igg
plot(zph_ba1_igg[1])
abline(h=ba1_igg_d28$coefficients[1], col=3, lwd=2, lty=2)

zph_ba1_iga <- cox.zph(ba1_iga_d28, transform = "identity")
zph_ba1_iga
plot(zph_ba1_iga[1])
abline(h=ba1_iga_d28$coefficients[1], col=3, lwd=2, lty=2)

wt_igg <- timecox(Surv(tstart, tstop, event) ~ antibody + Gender + age_cat 
                  + Booster + log2(cases), data = full_wt_igg, id = id, 
                  n.sim = 500)
wt_iga <- timecox(Surv(tstart, tstop, event) ~ antibody + Gender + age_cat 
                  + Booster + log(cases), data = full_wt_iga, id = id, 
                  n.sim = 500)
ba1_igg <- timecox(Surv(tstart, tstop, event) ~ antibody + Gender + age_cat 
                  + Booster + log(cases), data = full_ba1_igg, id = id, 
                  n.sim = 500)
ba1_iga <- timecox(Surv(tstart, tstop, event) ~ antibody + Gender + age_cat 
                  + Booster + log(cases), data = full_ba1_iga, id = id, 
                  n.sim = 500)
summary(wt_igg)
summary(wt_iga)
summary(ba1_igg)
summary(ba1_iga)
par(mfrow=c(2,3))
plot(wt_igg)
plot(wt_iga)
plot(ba1_igg)
plot(ba1_iga)

# Day 28;  25-50-25 KM -------------------------------------------------------------------
d28_cat4_wt_igg <- survfit(Surv(doi_new, Infect) ~ D28_antibody_level, data = d28_wt_igg_df)
wt_igg_km <- ggsurvplot(d28_cat4_wt_igg, risk.table = F, title = "WT IgG", 
           legend = "none", conf.int = F)

d28_cat4_wt_iga <- survfit(Surv(doi_new, Infect) ~ D28_antibody_level, data = d28_wt_iga_df)
wt_iga_km <- ggsurvplot(d28_cat4_wt_iga, risk.table = F, title = "WT IgA", 
           legend = "none", conf.int = F)

d28_cat4_ba1_igg <- survfit(Surv(doi_new, Infect) ~ D28_antibody_level, data = d28_ba1_igg_df)
ba1_igg_km <- ggsurvplot(d28_cat4_ba1_igg, risk.table = F, title = "BA.1 IgG", 
           legend = "none", conf.int = F)

d28_cat4_ba1_iga <- survfit(Surv(doi_new, Infect) ~ D28_antibody_level, data = d28_ba1_iga_df)
ba1_iga_km <- ggsurvplot(d28_cat4_ba1_iga, risk.table = F, title = "BA.1 IgA", 
                         legend = c(0.8,0.9), legend.title = "Category", 
                         legend.labs = c("Medium", "High", "Low"), conf.int = F)

# Plot all curves on same page
all_km <- list(a = wt_igg_km, b = ba1_iga_km)
arrange_ggsurvplots(all_km, print = TRUE, ncol = 2, nrow = 1)

# Survival for continuous antibody levels -------------------------------
wt_igg_cont <- coxph(Surv(tstart, tstop, event) ~ antibody 
                    + log2(cases) 
                    + Gender + age_cat + Booster, data = full_wt_igg) 
wt_iga_cont <- coxph(Surv(tstart, tstop, event) ~ antibody + 
                      log2(cases) +
                      Gender + age_cat + Booster, data = full_wt_iga) 
ba1_igg_cont <- coxph(Surv(tstart, tstop, event) ~ antibody + 
                       log2(cases) +
                       Gender + age_cat + Booster, data = full_ba1_igg) 
ba1_iga_cont <- coxph(Surv(tstart, tstop, event) ~ antibody + 
                       log2(cases) 
                      + Gender + age_cat + Booster, data = full_ba1_iga) 
summary(wt_igg_cont)
summary(wt_iga_cont)
summary(ba1_igg_cont)
summary(ba1_iga_cont) 

# Survival for categorical antibody levels -------------------------------
wt_igg_cat <- coxph(Surv(tstart, tstop, event) ~ cat4 
                    + log2(cases) 
                    + Gender + age_cat + Booster, data = full_wt_igg
                    #, maxit = 500,maxstep = 0.3
                    ) 
wt_iga_cat <- coxph(Surv(tstart, tstop, event) ~ cat4 + 
                      log2(cases) +
                      Gender + age_cat + Booster, data = full_wt_iga
                     #, maxit = 500, maxstep = 0.3
                     ) 
ba1_igg_cat <- coxph(Surv(tstart, tstop, event) ~ cat4 + 
                       log2(cases) +
                       Gender + age_cat + Booster, data = full_ba1_igg
                      #, maxit = 500, maxstep = 0.3
                      ) 

full_ba1_iga$cat4 <- 
  factor(full_ba1_iga$cat4, c("Low25%", "25-75%", "High25%"))
ba1_iga_cat <- coxph(Surv(tstart, tstop, event) ~ cat4 + 
                       log2(cases) +
                       Gender + age_cat + Booster, data = full_ba1_iga
                     #, maxit = 500, maxstep = 0.3
)

# Survival for categorical day 28 antibody levels-----------------------------------
d28_wt_igg <- coxph(Surv(doi_new, Infect) ~ D28_antibody_level + Gender + age_cat + Booster
                    + log2(cases)
                    , data = d28_wt_igg_df)
d28_wt_iga <- coxph(Surv(doi_new, Infect) ~ D28_antibody_level + Gender + age_cat + Booster
                    + log2(cases)
                    , data = d28_wt_iga_df)
d28_ba1_igg <- coxph(Surv(doi_new, Infect) ~ D28_antibody_level + Gender + age_cat + Booster
                     + log2(cases)
                     , data = d28_ba1_igg_df)
d28_ba1_iga <- coxph(Surv(doi_new, Infect) ~ D28_antibody_level + Gender + age_cat + Booster
                     + log2(cases)
                     , data = d28_ba1_iga_df)
summary(d28_wt_igg)
summary(d28_wt_iga)
summary(d28_ba1_igg)
summary(d28_ba1_iga)

# model comparison---------
d28_ab <- full_d28_ab %>% group_by(Code) %>% slice(1)

# Define time points to evaluate the ROC curves
time_points <- c(159, 249, 337)

# Time-dependent ROC for baseline model
roc_baseline <- timeROC(T = d28_ab$DoI, delta = d28_ab$inf, 
                        marker = predict(ba1_iga_d28), cause = 1, 
                        times = time_points, iid = F)

# Time-dependent ROC for time-varying model
roc_timevarying <- timeROC(T = full_ba1_iga$tstart, delta = full_ba1_iga$Infect, 
                           marker = predict(ba1_iga_cont), cause = 1, 
                           times = time_points, iid = F)

# Plot ROC curves for comparison
par(mfrow=c(3,2))
plot(roc_baseline, time = 159, col = "red", lwd = 2)
legend("bottomright", legend = "Baseline Model", col = "red", lwd = 2)
plot(roc_timevarying, time = 159, col = "blue", add = F, lwd = 2)
legend("topleft", legend = "Time-Varying Model", col = "blue", lwd = 2)

plot(roc_baseline, time = 249, col = "red", lwd = 2)
legend("bottomright", legend = "Baseline Model", col = "red", lwd = 2)
plot(roc_timevarying, time = 249, col = "blue", add = F, lwd = 2)
legend("topleft", legend = "Time-Varying Model", col = "blue", lwd = 2)

plot(roc_baseline, time = 337, col = "red", lwd = 2)
legend("bottomright", legend = "Baseline Model", col = "red", lwd = 2)
plot(roc_timevarying, time = 337, col = "blue", add = F, lwd = 2)
legend("topleft", legend = "Time-Varying Model", col = "blue", lwd = 2)

AIC(ba1_iga_d28)
AIC(ba1_iga_cont)
BIC(ba1_iga_d28)
BIC(ba1_iga_cont)
summary(ba1_iga_d28)
summary(ba1_iga_cont)

AIC(wt_igg_d28)
AIC(wt_igg_cont)
BIC(wt_igg_d28)
BIC(wt_igg_cont)
summary(wt_igg_d28)
summary(wt_igg_cont)