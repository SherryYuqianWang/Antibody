# Population fit ##############################################################
#generate samples from population parameter 
sample_ab <- function(pop, num, booster, infect_late,infect_un, age, gender){
  
  inds_mean <- which(rownames(pop) %in% c("k1_pop","k2_pop","k3_pop","k4_pop","A0_pop"))
  inds_sd <- which(rownames(pop) %in% c("omega_k1","omega_k2","omega_k3","omega_k4","omega_A0"))

  pars <- matrix(0, num+1, length(inds_mean)+1)
  
  for (i in 1:length(inds_mean)) {
    mean_par <- pop$value[inds_mean[i]]
    sd_par <- pop$value[inds_sd[i]]
    #sd_par <- pop$se_sa[inds_mean[i]] ####standard error of population parameter
    
    if (i==1 & (inds_mean[i+1]-inds_mean[i]==3)) {
      mean_infect_late <- pop$value[inds_mean[i]+1]
      mean_infect_un <- pop$value[inds_mean[i]+2]
      meanlog = log(mean_par)+infect_late*mean_infect_late+infect_un*mean_infect_un
      mean_log_k1 <- meanlog
      pars[1:num, i] <- exp(rnorm(num, mean = meanlog, sd=sd_par))
      pars[num+1, i] <- exp(meanlog)
    } 
    
    else if (i==1 & (inds_mean[i+1]-inds_mean[i]==2)) {
      mean_booster_pf <- pop$value[inds_mean[i]+1]
      meanlog = log(mean_par)+booster*mean_booster_pf
      mean_log_k1 <- meanlog
      pars[1:num, i] <- exp(rnorm(num, mean = meanlog, sd=sd_par))
      pars[num+1, i] <- exp(meanlog)
    } 
    
    else if (i==2) {
      mean_booster_pf <- pop$value[inds_mean[i]+1]
      meanlog = log(mean_par)+booster*mean_booster_pf
      pars[1:num, i] <- exp(rnorm(num, mean=meanlog, sd=sd_par))
      pars[num+1, i] <- exp(meanlog)
    }
    
    else if (i==3 & (inds_mean[i+1]-inds_mean[i]==4)) {
      mean_infect_late <- pop$value[inds_mean[i]+1]
      mean_booster_pf <- pop$value[inds_mean[i]+3]
      meanlog = log(mean_par)+infect_late*mean_infect_late+booster*mean_booster_pf
      pars[1:num, i] <- exp(rnorm(num, mean=meanlog, sd=sd_par))
      pars[num+1, i] <- exp(meanlog)
    }
    
    else if (i==3 & (inds_mean[i+1]-inds_mean[i]==5)) {
      mean_gender_male <- pop$value[inds_mean[i]+1]
      mean_infect_late <- pop$value[inds_mean[i]+2]
      mean_age <- pop$value[inds_mean[i]+4]
      meanlog = log(mean_par)+gender*mean_gender_male+infect_late*mean_infect_late+age*mean_age
      pars[1:num, i] <- exp(rnorm(num, mean=meanlog, sd=sd_par))
      pars[num+1, i] <- exp(meanlog)
    }
    
    else if (i==3 & (inds_mean[i+1]-inds_mean[i]==6)) {
      mean_gender_male <- pop$value[inds_mean[i]+1]
      mean_infect_late <- pop$value[inds_mean[i]+2]
      mean_booster_pf <- pop$value[inds_mean[i]+4]
      mean_age <- pop$value[inds_mean[i]+5]
      meanlog = log(mean_par)+gender*mean_gender_male+infect_late*mean_infect_late+booster*mean_booster_pf+age*mean_age
      pars[1:num, i] <- exp(rnorm(num, mean=meanlog, sd=sd_par))
      pars[num+1, i] <- exp(meanlog)
    }
    
    else if (i==4 & (rownames(pop)[inds_mean[i]+1] == "beta_k4__1st_Booster_vaccine_type_Pfizer")) {
      mean_booster_pf <- pop$value[inds_mean[i]+1]
      meanlog = log(mean_par)+booster*mean_booster_pf
      pars[1:num, i] <- exp(rnorm(num, mean=meanlog, sd=sd_par))
      pars[num+1, i] <- exp(meanlog)
    }
    
    else if (i==4 & (rownames(pop)[inds_mean[i]+1] == "beta_k4_Gender_Male")) {
      mean_gender_male <- pop$value[inds_mean[i]+1]
      meanlog = log(mean_par)+gender*mean_gender_male
      pars[1:num, i] <- exp(rnorm(num, mean=meanlog, sd=sd_par))
      pars[num+1, i] <- exp(meanlog)
    }
    
    else if (i==4 & (inds_mean[i+1]-inds_mean[i]==3)) {
      mean_infect_late <- pop$value[inds_mean[i]+1]
      meanlog = log(mean_par)+infect_late*mean_infect_late
      pars[1:num, i] <- exp(rnorm(num, mean=meanlog, sd=sd_par))
      pars[num+1, i] <- exp(meanlog)
    }
    
    else if (i==5 & (rownames(pop)[inds_mean[i]+1] == "beta_A0__1st_Booster_vaccine_type_Pfizer")) {
      mean_booster_pf <- pop$value[inds_mean[i]+1]
      meanlog = log(mean_par)+booster*mean_booster_pf
      mean_log_a0 <- meanlog
      pars[1:num, i] <- exp(rnorm(num, mean=meanlog, sd=sd_par))
      pars[num+1, i] <- exp(meanlog)
    }
    
    else { 
      mean_log_a0 <- log(mean_par)
      pars[1:num, i] <- exp(rnorm(num, mean=log(mean_par), sd=sd_par))
      pars[num+1, i] <- mean_par
    }
  }
  
  sd_k1 <- pop$value[inds_sd[1]]
  sd_a0 <- pop$value[inds_sd[5]]
  corr <- pop$value[which(row.names(pop)=="corr_k1_A0")]
  cov_k1_a0 <- corr*sd_k1*sd_a0
  sigma <- rbind(c(sd_k1^2, cov_k1_a0), c(cov_k1_a0, sd_a0^2))
  mu <- c(mean_log_k1, mean_log_a0)
  res <- mvrnorm(n=num, mu=mu, Sigma=sigma)
  pars[1:num, c(1, 5)] <- exp(res)
  return(pars)
}

#solve ode
ode_ab_pop<-function(pars,infect){
  
  k1 <- as.numeric(pars[1])
  k2 <- as.numeric(pars[2])
  k3 <- as.numeric(pars[3])
  k4 <- as.numeric(pars[4])
  A0 <- as.numeric(pars[5])
  
  times<-c(seq(Tmin,Tmax,step_size))
  infection <- data.frame(times = times, inf_0 = rep(0, length(times)), inf_1 = rep(0, length(times)),inf_2 = rep(0, length(times)))
  infection <- infection %>% mutate(inf_0 =
                                     case_when(times <= 21 ~ k1, 
                                               times > 21 & times <= 112 ~ -k2,
                                               times > 112 & times <= 133 ~ k3,
                                               times > 133 ~ -k4), 
                                   inf_1 = 
                                     case_when(times <= 21 ~ k1, 
                                               times > 21 & times <= 264 ~ -k2,
                                               times > 264 & times <= 285 ~ k3,
                                               times > 285 ~ -k4),
                                   inf_2 = 
                                     case_when(times <= 21 ~ k1, 
                                               times > 21 ~ -k2))
  
  k_t <- approxfun(infection$times, infection[,infect+2], rule = 2)
  
  derivs<-function(times,y,pars,k_t){
    with(as.list(c(pars,y)),{
      k <- k_t(times)   # <---- here
      dA <- k*A
      return(list(c(dA)))
    })
  }
  y<-c(A=A0)
  
  
  out<-ode(y=y,parms=pars,times=times,func=derivs, k_t=k_t)
  as.data.frame(out)
}

###simulation to plot the confidence interval
simulate_AB<-function(group_list){
  
  booster <- group_list$booster
  infect_late <- group_list$infect_late
  infect_un <- group_list$infect_un
  infect <- group_list$infect
  age <- group_list$age
  gender <- group_list$gender
  
  pars <- sample_ab(pop, num, booster, infect_late, infect_un, age, gender)
  
  ind_k1 <- as.numeric(pars[,1])
  ind_k2 <- as.numeric(pars[,2])
  ind_k3 <- as.numeric(pars[,3])
  ind_k4 <- as.numeric(pars[,4])
  ind_A0 <- as.numeric(pars[,5])
  
  total_cov <- as.data.frame(cbind(booster,infect_late,infect_un,infect,age,gender,ind_k1,ind_k2,ind_k3,ind_k4,ind_A0))
  
  AB <- matrix(NA,nrow=length(seq(Tmin,Tmax,step_size)),ncol=nrow(total_cov))
  
  for(i in 1:nrow(total_cov)){
    
    par <- c(k1=total_cov$ind_k1[i],
             k2=total_cov$ind_k2[i],
             k3=total_cov$ind_k3[i],
             k4=total_cov$ind_k4[i],
             A0=total_cov$ind_A0[i])
    out <- ode_ab_pop(par,infect)
    AB[, i] <- out$A
  }
  
  ######calculate mean viral load
  mean_ab  <- apply(AB,1,function(x){quantile(x, 0.5, na.rm=T)})
  low_ab  <- apply(AB,1,function(x){quantile(x, 0.05, na.rm=T)})
  q1_ab  <- apply(AB,1,function(x){quantile(x, 0.25, na.rm=T)})
  q3_ab  <- apply(AB,1,function(x){quantile(x, 0.75, na.rm=T)})
  high_ab  <- apply(AB,1,function(x){quantile(x, 0.95, na.rm=T)})
  best_ab <- AB[,nrow(total_cov)]
  
  Fit <- cbind(seq(Tmin,Tmax,step_size),mean_ab,low_ab,q1_ab,q3_ab,high_ab,best_ab,booster,infect,age,gender)
  Fit <- data.frame(Fit)
  Fit$q3_ab[Fit$q3_ab > 70] <- 70
  Fit$high_ab[Fit$high_ab > 70] <- 70
  
  colnames(Fit) <- c("time","mean_ab","low_ab","q1_ab","q3_ab","high_ab","best_ab","Booster","Infect","Age","Gender")
  return(Fit)
}



# Individual fit ###############################################################
#solve ode (individual fit with k1,k2,k3,k4) 
ode_ab_ind<-function(pars,infect,doi){
  
  k1 <- as.numeric(pars[1])
  k2 <- as.numeric(pars[2])
  k3 <- as.numeric(pars[3])
  k4 <- as.numeric(pars[4])
  A0 <- as.numeric(pars[5])
  
  times<-c(seq(Tmin,Tmax,step_size))
  infection <- data.frame(times = times, inf_0 = rep(0, length(times)), inf_1 = rep(0, length(times)),inf_2 = rep(0, length(times)))
  infection <- infection %>% mutate(inf_0 =
                                      case_when(times <= 21 ~ k1, 
                                                times > 21 & times <= doi ~ -k2,
                                                times > doi & times <= doi+21 ~ k3,
                                                times > doi+21 ~ -k4), 
                                    inf_1 = 
                                      case_when(times <= 21 ~ k1, 
                                                times > 21 & times <= doi ~ -k2,
                                                times > doi & times <= doi+21 ~ k3,
                                                times > doi+21 ~ -k4),
                                    inf_2 = 
                                      case_when(times <= 21 ~ k1, 
                                                times > 21 ~ -k2))
  
  k_t <- approxfun(infection$times, infection[,infect+2], rule = 2)
  
  derivs<-function(times,y,pars,k_t){
    with(as.list(c(pars,y)),{
      k <- k_t(times)   # <---- here
      dA <- k*A
      return(list(c(dA)))
    })
  }
  y<-c(A=A0)
  
  
  out<-ode(y=y,parms=pars,times=times,func=derivs, k_t=k_t)
  as.data.frame(out)
}

#solve ode (individual fit with k1,k2) 
ode_ab_ind_cluster<-function(pars,infect){
  
  k1 <- as.numeric(pars[1])
  k2 <- as.numeric(pars[2])
  #k3 <- as.numeric(pars[3])
  #k4 <- as.numeric(pars[4])
  A0 <- as.numeric(pars[5])
  
  times<-c(seq(Tmin,Tmax,step_size))
  infection <- data.frame(times = times, inf_0 = rep(0, length(times)), inf_1 = rep(0, length(times)),inf_2 = rep(0, length(times)))
  infection <- infection %>% mutate(inf_0 =
                                      case_when(times <= 21 ~ k1,
                                                times > 21 ~ -k2), 
                                    inf_1 = 
                                      case_when(times <= 21 ~ k1, 
                                                times > 21 ~ -k2),
                                    inf_2 = 
                                      case_when(times <= 21 ~ k1, 
                                                times > 21 ~ -k2))
  
  k_t <- approxfun(infection$times, infection[,infect+2], rule = 2)
  
  derivs<-function(times,y,pars,k_t){
    with(as.list(c(pars,y)),{
      k <- k_t(times)   # <---- here
      dA <- k*A
      return(list(c(dA)))
    })
  }
  y<-c(A=A0)
  
  
  out<-ode(y=y,parms=pars,times=times,func=derivs, k_t=k_t)
  as.data.frame(out)
}

## Table ####
#generate dataset with k1,k2,k3,k4
ind_ab_long <- function(Estimated, original_cov, sg_covid){
  
  Est <- Estimated %>% mutate(Infect =
                                case_when(Infection == "Early" ~ 0,
                                          Infection == "Late" ~ 1,
                                          Infection == "Uninfected" ~ 2),
                              Booster = 
                                case_when(X_1st_Booster_vaccine_type == "Moderna" ~ 0,
                                          X_1st_Booster_vaccine_type == "Pfizer" ~ 1))
  
  Est <- merge(Est, original_cov, by.x = "id", by.y = "Code")
  Est <- Est %>%
    add_column(booster_peak = NA,
               infection_peak = NA,
               infect_antibody = NA,
               censor = NA,
               infection_fold = NA)
  
  AB <- matrix(NA,nrow=nrow(Est),ncol=length(seq(Tmin+22,Tmax,step_size)))
  
  colname <- paste("antibody",c(1:(Tmax-21)))
  colnames(AB) <- colname
  
  case_weekly <- matrix(NA,nrow=nrow(Est),ncol=length(seq(Tmin+22,Tmax,step_size)))
  
  colname <- paste("cases",c(1:(Tmax-21)))
  colnames(case_weekly) <- colname
  
  
  df_antibody<-cbind(Est, AB, case_weekly)
  
  for(i in 1:nrow(Est)){
    infect <- Est$Infect[i]
    doi <- Est$DoI[i]
    b <- which(sg_covid[,1] == Est$booster_trans[i])
    #adjust <-Est$adjust[i]
    #booster <- Est$Booster[i]
    pars <- c(k1=Est$k1_SAEM[i],
              k2=Est$k2_SAEM[i],
              k3=Est$k3_SAEM[i],
              k4=Est$k4_SAEM[i],
              A0=Est$A0_SAEM[i])
    fitted <- ode_ab_ind(pars,infect,doi)
    d1 <- data.frame(x=(times),y=(fitted$A))
    Est$booster_peak[i] <- d1$y[23]  #peak at day 22 (row 23)
    Est$infection_peak[i] <- d1$y[doi+23]
    Est$infection_fold[i] <- d1$y[doi+23]/d1$y[doi+1]
    
    if (infect == 2){
      #Est$infect_antibody[i] <- d1$y[361]
      Est$infect_antibody[i] <- 0
      Est$censor[i] <- 1
    }
    else {
      Est$infect_antibody[i] <- d1$y[doi+1]
      Est$censor[i] <- 0
    }
    #df_antibody[i,]<-cbind(Est[i,],t(d1$y[(adjust+1):(adjust+339)])) ####adjust to same day 0
    df_antibody[i,]<-cbind(Est[i,],t(d1$y[23:(Tmax+1)]),t(sg_covid$Singapore[(b+22):(b+Tmax)]))  ####different day 0
    
  }
  
  antibody_long <- subset(df_antibody, select = -c(2:21)) %>% 
    pivot_longer(cols = -c(1:15), names_to = c(".value", "time"), names_sep = " ") %>%
    mutate_at(
      .vars = "time",
      .funs = as.numeric
    )
  
  colnames(antibody_long)[10] <- "day28_ab"
  
  ind_AB <- antibody_long %>% 
    mutate(doi_new = ifelse(is.na(DoI), 360, DoI-21))%>%
    mutate(antibody_new = ifelse(time<=doi_new, antibody, as.integer(NA)))%>%
    mutate(antibody_cat4 = ntile(antibody_new, 4),
           day28_cat4 = ntile(day28_ab,4)) %>%
    mutate(antibody_cat4 = case_when(antibody_cat4 == 1 ~ "Low25%",
                                     antibody_cat4 == 2 | antibody_cat4 == 3 ~ "25-75%",
                                     antibody_cat4 == 4 ~ "High25%"),
           day28_cat4 = case_when(day28_cat4 == 1 ~ "Low25%",
                                  day28_cat4 == 2 | day28_cat4 == 3 ~ "25-75%",
                                  day28_cat4 == 4 ~ "High25%")) %>%
    mutate(doi_new = ifelse(is.na(DoI), as.integer(NA), DoI-21))
  
  return(ind_AB)
}

#generate dataset with k1,k2
ind_ab_long_cluster <- function(Estimated, original_cov){
  
  Est <- Estimated[[4]]%>% mutate(Infect =
                               case_when(Infection == "Early" ~ 0,
                                         Infection == "Late" ~ 1,
                                         Infection == "Uninfected" ~ 2),
                             Booster = 
                               case_when(X_1st_Booster_vaccine_type == "Moderna" ~ 0,
                                         X_1st_Booster_vaccine_type == "Pfizer" ~ 1),
                             age_cat = 
                               case_when(age_cat == "<60" ~ "Age<60",
                                         age_cat == ">=60" ~ "Age>=60"))
  
  Est <- merge(Est, original_cov, by.x = "id", by.y = "Code")
  Est <- Est %>%
    add_column(booster_peak = NA,
               infection_peak = NA,
               infect_antibody = NA,
               censor = NA,
               infection_fold = NA,
               PT = NA,
               protection = NA,
               protection_i = NA)
  
  AB <- matrix(NA,nrow=nrow(Est),ncol=length(seq((Tmin+22),Tmax,step_size)))
  
  colname <- paste("antibody",c(22:Tmax))
  colnames(AB) <- colname
  
  df_antibody<-cbind(Est, AB)
  df_antibody_i<-cbind(Est, AB)
  
  for(i in 1:nrow(Est)){
    
    infect <- Est$Infect[i]
    doi <- Est$DoI[i]
    pars <- c(k1=Est$k1_mode[i],
              k2=Est$k2_mode[i],
              k3=Est$k3_mode[i],
              k4=Est$k4_mode[i],
              A0=Est$A0_mode[i])
    fitted_i <- ode_ab_ind(pars,infect,doi)
    fitted <- ode_ab_ind_cluster(pars,infect)
    
    d1 <- data.frame(x=(times),y=(fitted$A))
    
    if (is.na(doi)){
      d2 <- data.frame(x=(times),z=(fitted_i$A[1:(Tmax+1)]))
    } else {
      d2 <- data.frame(x=(times),y=(fitted$A),z=(fitted_i$A[(doi+1):(doi+Tmax+1)]))
    }    
    
    Est$booster_peak[i] <- d1$y[23]  #peak at day 22 (row 23)
    
    
    if (infect == 2){
      Est$infection_peak[i] <- NA
      Est$infection_fold[i] <- NA
      Est$infect_antibody[i] <- 0
      Est$censor[i] <- 1
    } else {
      Est$infection_peak[i] <- d2$z[23]
      Est$infection_fold[i] <- d2$z[23]/d2$z[1]
      Est$infect_antibody[i] <- d2$z[1]
      Est$censor[i] <- 0
    }
    
    df_antibody[i,]<-cbind(Est[i,],t(d1$y[23:(Tmax+1)]))
    df_antibody_i[i,]<-cbind(Est[i,],t(d2$z[23:(Tmax+1)]))
  }
  
  PT <- 10 ###80% protection in 3 months (only for BA1.IgA)
  
  for(i in 1:nrow(df_antibody)){
    df_antibody$PT[i] <- PT
    doi <- df_antibody$DoI[i]
    
    if (df_antibody[i,(Tmax+38-22)]>PT){ #column 38 is day23
      df_antibody$protection[i] <- 361
    } else {
      df_antibody$protection[i] <- min(which(df_antibody[i,38:(Tmax+38-22)] <= PT))+21
    }
  }
  
  antibody_long <- subset(df_antibody, select = -c(2:21)) %>% 
    pivot_longer(cols = -c(1:18), names_to = c(".value", "time"), names_sep = " ") %>%
    mutate_at(
      .vars = "time",
      .funs = as.numeric
    )
  
  colnames(antibody_long)[10] <- "day28_ab"
  
  ind_AB <- antibody_long %>% 
    mutate(antibody_new = ifelse(time<=DoI | is.na(DoI), antibody, as.integer(NA)))%>%
    mutate(antibody_cat4 = ntile(antibody_new, 4),
           day28_cat4 = ntile(day28_ab,4)) %>%
    mutate(antibody_cat4 = case_when(antibody_cat4 == 1 ~ "Low25%",
                                     antibody_cat4 == 2 | antibody_cat4 == 3 ~ "25-75%",
                                     antibody_cat4 == 4 ~ "High25%"),
           day28_cat4 = case_when(day28_cat4 == 1 ~ "Low25%",
                                  day28_cat4 == 2 | day28_cat4 == 3 ~ "25-75%",
                                  day28_cat4 == 4 ~ "High25%"))
  
  return(ind_AB)
}

## Figure ####
ind_fit_plt<-function(Estimated,Simulated,original_cov){
  
  Est <- Estimated %>% mutate(Infect =
                                case_when(Infection == "Early" ~ 0,
                                          Infection == "Late" ~ 1,
                                          Infection == "Uninfected" ~ 2))
  
  Est <- merge(Est, original_cov, by.x = "id", by.y = "Code")
  Fit <- list()
  for(i in 1:nrow(Est)){
    infect <- Est$Infect[i]
    doi <- Est$DoI[i]
    pars <- c(k1=Est$k1_SAEM[i],
              k2=Est$k2_SAEM[i],
              k3=Est$k3_SAEM[i],
              k4=Est$k4_SAEM[i],
              A0=Est$A0_SAEM[i])
    fitted <- ode_ab_ind(pars,infect,doi)
    d1 <- data.frame(Day=(times),y=(fitted$A))
    Code <- Est$id[i]
    gender <- Est$Gender[i]
    booster <- Est$X_1st_Booster_vaccine_type[i]
    age_cat <- Est$age_cat[i]
    
    S <- 10 ###10 repeats
    P <- matrix(NA,nrow=length(times),ncol=S)
    
    
    for(j in 1:S){
      pars <- c(k1=Simulated$k1[j+S*(i-1)],k2=Simulated$k2[j+S*(i-1)],k3=Simulated$k3[j+S*(i-1)],k4=Simulated$k4[j+S*(i-1)],A0=Simulated$A0[j+S*(i-1)])
      out  <- ode_ab_ind(pars,infect,doi)
      P[,j] <- out$A
    }
    
    Min  <- apply(P,1,function(x){quantile(x,0.025)})
    Max  <- apply(P,1,function(x){quantile(x,0.975)})
    
    fit <- cbind(d1,Min,Max,Code,gender,booster,age_cat,infect,doi)
    fit$Max[fit$Max> 100] <- 100
    fit$Min[fit$Min> 100] <- 100
    
    Fit[[i]] <- data.frame(fit)
  }
  
  ind_fit_unlist <- map_df(Fit, ~as.data.frame(.x))
  ind_fit <- merge(ind_fit_unlist, original[,c(1,10,11,12,13,14)], by=c("Code", "Day"),all=TRUE) 
  ind_fit$Code <- as.factor(ind_fit$Code)
  colnames(ind_fit)[11:14] <- c("WT IgG","BA1 IgG","WT IgA","BA1 IgA")
  return(ind_fit)
}



# Regression ###################################################################
#generate regression table
regression_table<-function(ind_AB){
  
  df_reg <- ind_AB %>% group_by(id) %>% slice(7)
  colnames(df_reg)[4] <- "booster"
  
  explanatory = c("Gender", "booster", "age_cat","Infection")
  dependent_all = c("day28_ab","booster_peak","infection_peak","protection")
  
  t_reg <- list()
  
  for (i in 1:3){
    dependent <- dependent_all[i]
    df_reg %>%
      ## Crosstable
      summary_factorlist(dependent, explanatory, fit_id=TRUE)  %>% 
      
      ## Add univariable
      ff_merge(
        lmuni(df_reg, dependent, explanatory) %>%
          fit2df(estimate_suffix=" (univariable)",last_merge = TRUE)
      ) %>% 
      
      ## Add multivariable
      ff_merge(
        lmmulti(df_reg, dependent, explanatory) %>%
          fit2df(estimate_suffix="",last_merge = TRUE)
      ) %>% 
      
      dependent_label(df_reg, dependent) %>%
      ff_remove_p() -> t_reg[[i]]
    
    t_reg[[i]] <- t_reg[[i]]#[,c(1,2,3,7,8,9)]
  }
  
  t_all <- reduce(.x = t_reg, merge, by = c('fit_id'), all = T)
  return(t_all)
}



# Survival analysis ############################################################
###
ode_ab_surv<-function(pars,tstop_ind){
  
  A0 <- pars[1]
  k2 <- pars[2]
  
  times<-c(seq(Tmin,tstop_ind,step_size))
  infection <- data.frame(times = times, inf_0 = rep(-k2, length(times)))
  #infection <- infection %>% mutate(inf_0 =
  #                                   case_when(times <= 21 ~ k1, 
                                                #times > 21 & times <= boost ~ -k2,
                                                #times > boost & times <= boost+21 ~ k3,
                                                #times > boost+21 ~ -k4))
  
  k_t <- approxfun(infection$times, infection[,2], rule = 2)
  
  derivs<-function(times,y,pars,k_t){
    with(as.list(c(pars,y)),{
      k <- k_t(times)   # <---- here
      dA <- k*A
      return(list(c(dA)))
    })
  }
  y<-c(A=A0)
  
  
  out<-ode(y=y,parms=pars,times=times,func=derivs, k_t=k_t)
  as.data.frame(out)
}


###
ode_ab_boost<-function(pars,boost_interval,tstop_ind){
  
  A0 <- pars[1] #antibody peak
  k1 <- pars[2]
  k2 <- pars[3]
  k3 <- pars[4]
  
  times<-c(seq(Tmin,tstop_ind,step_size))
  infection <- data.frame(times = times, inf_0 = c(rep(k1, 21),rep(-k2, length(times)-21)))
  
  #if (boost_interval==3){
  #  infection <- infection %>% mutate(inf_0 = ifelse((times>69 & times <90) | (times>159 & times <180) | (times>249 & times <270), k3,-k2))
  #} else 
    if (boost_interval==0){
    infection <- infection %>% mutate(inf_0 = 
                                        case_when(times <= 21 ~ k1, 
                                                  times > 21 & times <= 180 ~ -k2,
                                                  times > 180 & times <= 201 ~ k3))
  } else if (boost_interval==1){
    infection <- infection %>% mutate(inf_0 = 
                                        case_when(times <= 21 ~ k1, 
                                                  times > 21 & times <= 270 ~ -k2,
                                                  times > 270 & times <= 291 ~ k3))
  } 

  k_t <- approxfun(infection$times, infection[,2], rule = 2)
  
  derivs<-function(times,y,pars,k_t){
    with(as.list(c(pars,y)),{
      k <- k_t(times)   # <---- here
      dA <- k*A
      return(list(c(dA)))
    })
  }
  y<-c(A=A0)
  
  
  out<-ode(y=y,parms=pars,times=times,func=derivs, k_t=k_t)
  as.data.frame(out)
}
