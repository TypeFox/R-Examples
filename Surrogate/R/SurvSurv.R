SurvSurv <- function(Dataset, Surr, SurrCens, True, TrueCens, Treat, 
Trial.ID, Weighted=TRUE, Alpha=.05){ #, Number.Bootstraps=500, Seed=sample(1:1000, size=1)

  
  Surr <- Dataset[,paste(substitute(Surr))]
  SurrCens <- Dataset[,paste(substitute(SurrCens))] 
  True <- Dataset[,paste(substitute(True))]
  TrueCens <- Dataset[,paste(substitute(TrueCens))]
  Treat <- Dataset[,paste(substitute(Treat))]
  Trial.ID <- Dataset[,paste(substitute(Trial.ID))]
  All_data <- data.frame(cbind(Surr, SurrCens, True, TrueCens, Treat, Trial.ID))
  names(All_data) <- c("Surr", "SurrCens", "True", "TrueCens", "Treat", "Trial.ID")

  #R2h ind
  surv_data_Surr <- Surv(time=All_data$Surr, time2=All_data$SurrCens)
  surv_data_True <- Surv(time=All_data$True, time2=All_data$TrueCens)
  
  options(warn = -1)
  Model.0 <- 
    (coxph(surv_data_True~-1 + as.factor(All_data$Trial.ID) + 
             as.factor(All_data$Trial.ID):All_data$Treat, ties = "breslow"))
  Num.Events <- Model.0$nevent
  Model.1 <- 
    (coxph(surv_data_True~-1 + as.factor(All_data$Trial.ID) + 
             as.factor(All_data$Trial.ID):All_data$Treat +
             as.factor(All_data$Trial.ID):All_data$Surr, ties = "breslow"))
  L0 <- -2 * logLik(Model.0)[1]    
  L1 <- -2 * logLik(Model.1)[1]   
  
  if (logLik(Model.0)[1] > logLik(Model.1)[1]){
  warning("the likelihood of model True~Treat is higher than model True~Treat+Surr. Examine the fit of the models. R2_h.ind is unreliable. ")
  }
  
  # R2h.ind (single-cluster based) and CI
  g2 <- -(L1-L0)
  R2h.ind.value <- 1 - exp(-g2/length(All_data$Treat))
  options(warn=0)
  k1 <- qchisq(Alpha, 1, g2)
  d1 <- qchisq((1-Alpha), 1, g2)
  R2h.ind.lb <- max(0, 1-exp(-k1/length(All_data$Treat))) 
  R2h.ind.ub <- min(1, 1-exp(-d1/length(All_data$Treat)))
  R2h.ind <- data.frame(cbind(R2h.ind.value, R2h.ind.lb, R2h.ind.ub))   # OUT
  colnames(R2h.ind) <- c("R2h.ind", "CI lower limit", "CI upper limit")
  rownames(R2h.ind) <- c(" ") 
  
  # R2h.ind (single-cluster based) and CI
  g2 <- -(L1-L0)
  R2h.ind.value <- 1 - exp(-g2/Num.Events)
  options(warn=0)
  k1 <- qchisq(Alpha, 1, g2)
  d1 <- qchisq((1-Alpha), 1, g2)
  R2h.ind.lb <- max(0, 1-exp(-k1/Num.Events)) 
  R2h.ind.ub <- min(1, 1-exp(-d1/Num.Events))
  R2h.ind.QF <- data.frame(cbind(R2h.ind.value, R2h.ind.lb, R2h.ind.ub))   # OUT
  colnames(R2h.ind.QF) <- c("R2h.ind.QF", "CI lower limit", "CI upper limit")
  rownames(R2h.ind.QF) <- c(" ") 
  
  
  # R2_h indiv per trial
  R2h.per.trial.lb.vec <- R2h.per.trial.Name.vec <- R2h.per.trial.ub.vec <- R2h.per.trial.hier.vec <- NULL
  trialnames <- unique(All_data$Trial.ID)
  
  options(warn=-1)
  for (i in 1:length(unique(All_data$Trial.ID))) {
    Data.temp <- subset(All_data, subset=(Trial.ID==trialnames[i]))  # inf likelihoods als <4 observaties voor Model.1
    if (nrow(Data.temp) >= 4) {    
      if (length(unique(Data.temp$True))!=1){          # do not consider clusters with all same T (inf L)
        if (length(unique(Data.temp$Treat))==2){        # not consider clusters with all same Treat
          R2h.per.trial.Name <- trialnames[i]
          # R2h.per.trial 
          surv_data_Surr <- Surv(time=Data.temp$Surr, time2=Data.temp$SurrCens)
          surv_data_True <- Surv(time=Data.temp$True, time2=Data.temp$TrueCens)
          
          options(warn = -1)
          Model.0 <- 
            (coxph(surv_data_True~ 1 + Data.temp$Treat, ties = "breslow"))
          Num.Events <- Model.0$nevent
          Model.1 <- 
            (coxph(surv_data_True~ 1 + Data.temp$Treat + Data.temp$Surr, ties = "breslow"))
          L0 <- -2 * logLik(Model.0)[1]    
          L1 <- -2 * logLik(Model.1)[1]   
          
          if (logLik(Model.0)[1] > logLik(Model.1)[1]){
            warning("the likelihood of model True~Treat is higher than model True~Treat+Surr. Examine the fit of the models. R2_h.ind is unreliable. ")
          }
          
          # R2h.ind (single-cluster based) and CI
          g2.hier <- -(L1-L0)
            
          if (g2.hier != 0){
            R2h.per.trial.hier <- 1 - exp(-g2.hier/Num.Events)
            k1 <- qchisq(Alpha, 1, g2.hier)
            d1 <- qchisq((1-Alpha), 1, g2.hier)
            R2h.per.trial.lb <- max(0, 1-exp(-k1/Num.Events)) 
            R2h.per.trial.ub <- min(1, 1-exp(-d1/Num.Events))
            
            R2h.per.trial.hier.vec <- cbind(R2h.per.trial.hier.vec, R2h.per.trial.hier)
            R2h.per.trial.lb.vec <- cbind(R2h.per.trial.lb.vec, R2h.per.trial.lb)
            R2h.per.trial.ub.vec <- cbind(R2h.per.trial.ub.vec, R2h.per.trial.ub)
            R2h.per.trial.Name.vec <- cbind(R2h.per.trial.Name.vec, R2h.per.trial.Name)
          }}}}}
  options(warn=0)
  
  R2h.By.Trial.QF <- cbind(t(R2h.per.trial.Name.vec), t(R2h.per.trial.hier.vec), t(R2h.per.trial.lb.vec), 
                        t(R2h.per.trial.ub.vec))
  colnames(R2h.By.Trial.QF) <- c("TrialID", "R2h.ind", "R2h_low", "R2h_up")
  rownames(R2h.By.Trial.QF) <- NULL 
  R2h.By.Trial.QF <- data.frame(na.exclude(R2h.By.Trial.QF))      
  
  
  
  
  # Bootstrap
#  num <- dim(All_data)[1]
#  R2h.ind.value.boot.all <- NULL
#  for (i in 1: Number.Bootstraps){
#    obs <- c(1:num)
#    set.seed(Seed)
#    Indicator <- sample(obs, num, replace=TRUE)  #sample with replacement from 1...Ntotal
#    Sample.Boot <- data.frame(All_data[Indicator,])  
#    Seed <- Seed + 1
    
#    Data_boot_all <- NULL
    
    # remove trials that do not have >=3 observations in each treatment arm due to estimability constraints 
#    for (k in 1: length(unique(Sample.Boot$Trial.ID))){
#      Data_hier <- 
#        Sample.Boot[Sample.Boot$Trial.ID == unique(Sample.Boot$Trial.ID)[k],]
 
#    Data_hier$Treat[Data_hier$Treat==-1] <- 0
    
#    if ((dim(Data_hier[Data_hier$Treat==0,])[1]>=3) & (dim(Data_hier[Data_hier$Treat==1,])[1])>=3){
#    Data_boot_all <- rbind(Data_boot_all, Data_hier)}
#    }
    
#    surv_data_Surr <- Surv(time=Data_boot_all$Surr, time2=Data_boot_all$SurrCens)
#    surv_data_True <- Surv(time=Data_boot_all$True, time2=Data_boot_all$TrueCens)
    
#    Model.0.boot <- 
#      (coxph(surv_data_True~-1 + as.factor(Data_boot_all$Trial.ID) + 
#               as.factor(Data_boot_all$Trial.ID):Data_boot_all$Treat, ties = "breslow"))
#    Num.Events.boot <- Model.0.boot$nevent
#    Model.1.boot <- 
#      (coxph(surv_data_True~-1 + as.factor(Data_boot_all$Trial.ID) + 
#               as.factor(Data_boot_all$Trial.ID):Data_boot_all$Treat +
#               as.factor(Data_boot_all$Trial.ID):Data_boot_all$Surr, ties = "breslow"))
    
#    L0.boot <- -2 * logLik(Model.0.boot)[1]    
#    L1.boot <- -2 * logLik(Model.1.boot)[1]   
    
#    if (logLik(Model.0.boot)[1] <= logLik(Model.1.boot)[1]){
#    g2.boot <- -(L1.boot-L0.boot)
#    R2h.ind.value.boot <- 1 - exp(-g2.boot/Num.Events.boot)
#    R2h.ind.value.boot.all <- cbind(R2h.ind.value.boot.all, R2h.ind.value.boot)
#    }
#  }
# nonparametric bootstrap normal confidence interval
#  Mean.Boot.CI <- mean(R2h.ind.value.boot.all)
#  Var.Boot.CI <- var(Boot.CI)
#  Sort.CI <- sort(Boot.CI)
#  Boot.CI.R2.lb <- max(0, R2h.ind.c.b.val + qnorm(Alpha/2)*sqrt(Var.Boot.CI))    
#  Boot.CI.R2.ub <- min(1, R2h.ind.c.b.val - qnorm(Alpha/2)*sqrt(Var.Boot.CI))
#  R2h.cluster.based <- data.frame(cbind(R2h.ind.c.b.val, sqrt(Var.Boot.CI), Boot.CI.R2.lb, Boot.CI.R2.ub))
#  colnames(R2h.cluster.based) <- c("R2h", "Standard Error", "CI lower limit", "CI upper limit")
#  rownames(R2h.cluster.based) <- c(" ")
  

  # Trial
  Results.Stage.1 <- Data.Analyze <- NULL
  for (i in 1: length(unique(Trial.ID))){
    Data_hier <- All_data[All_data$Trial.ID == unique(Trial.ID)[i],]
    names(Data_hier) <- c("Surr", "SurrCens", "True", "TrueCens", "Treat", "Trial.ID")
    Trial.Size <- dim(Data_hier)[1]
    Trial.Name <- unique(Trial.ID)[i]
    Data_hier$Treat[Data_hier$Treat==-1] <- 0 
    
    if ((dim(Data_hier[Data_hier$Treat==0,])[1]<3) & ((dim(Data_hier[Data_hier$Treat==1,])[1])<3)){
     cat("\nNote. The trial with ID ", unique(Trial.ID)[i], " did not have >=3 observations in each treatment arm and \nwas excluded from the trial-level analyses (estimation of R2_ht) due to estimability constraints. \n", sep="")}    
    
    if ((dim(Data_hier[Data_hier$Treat==0,])[1]>=3) & ((dim(Data_hier[Data_hier$Treat==1,])[1])>=3)){
      surv_data_Surr <- Surv(time=Data_hier$Surr, time2=Data_hier$SurrCens)
      LogHazard_Surr <- as.numeric(coxph(surv_data_Surr~Data_hier$Treat, ties = "breslow")[1])
      surv_data_True <- Surv(time=Data_hier$True, time2=Data_hier$TrueCens)
      LogHazard_True <- as.numeric(coxph(surv_data_True~Data_hier$Treat, ties = "breslow")[1])

      Results_here <- cbind(Trial.Name, LogHazard_Surr, LogHazard_True, Trial.Size)
      Results.Stage.1 <- rbind(Results.Stage.1, Results_here)
      Data.Analyze <- rbind(Data.Analyze, Data_hier)
    }
  }
  Results.Stage.1 <- data.frame(Results.Stage.1)
  
  # stage 2
  if (Weighted==FALSE){
      Results.Stage.2 <- lm(Results.Stage.1$LogHazard_True ~ Results.Stage.1$LogHazard_Surr)
    }
  
  if (Weighted==TRUE){
      Results.Stage.2 <- lm(Results.Stage.1$LogHazard_True ~ Results.Stage.1$LogHazard_Surr, weights=Results.Stage.1$Trial.Size) 
  }
  
  # Trial R2
  N.trial <- length(unique(Results.Stage.1$Trial.Name)) 
  Trial.R2.value <- as.numeric(summary(Results.Stage.2)[c("r.squared")])
  Trial.R2.sd <- sqrt((4*Trial.R2.value*(1-Trial.R2.value)^2)/(N.trial-3))
  Trial.R2.lb <- max(0, Trial.R2.value + qnorm(Alpha/2) *(Trial.R2.sd))
  Trial.R2.ub <- min(1, Trial.R2.value + qnorm(1-Alpha/2)*(Trial.R2.sd))
  Trial.R2 <- data.frame(cbind(Trial.R2.value, Trial.R2.sd, Trial.R2.lb, Trial.R2.ub))
  colnames(Trial.R2) <- c("R2_trial", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(Trial.R2) <- c(" ") 
  
  # Trial R
  Trial.R.value <- sqrt(as.numeric(summary(Results.Stage.2)[c("r.squared")]))
  Z <- .5*log((1+Trial.R.value)/(1-Trial.R.value)) 
  Trial.R.lb <- max(0, (exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))-1)/(exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))+1))
  Trial.R.ub <- min(1, (exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))-1)/(exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))+1))
  Trial.R.sd <- sqrt((1-Trial.R.value**2)/(N.trial-2))
  Trial.R <- data.frame(cbind(Trial.R.value, Trial.R.sd , Trial.R.lb, Trial.R.ub))
  colnames(Trial.R) <- c("R_trial", "Standard Error", "CI lower limit", "CI upper limit")
  row.names(Trial.R) <- c(" ") 
  

  fit <- 
    list(Results.Stage.1=Results.Stage.1, 
         Results.Stage.2=Results.Stage.2, R2.trial=Trial.R2, R2.hind=R2h.ind, 
         R2h.ind.QF=R2h.ind.QF, R2.hInd.By.Trial.QF=R2h.By.Trial.QF, 
         Call=match.call())
  
  class(fit) <- "SurvSurv"
  fit
  
}


