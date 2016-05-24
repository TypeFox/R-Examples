TwoStageSurvSurv <- function(Dataset, Surr, SurrCens, True, TrueCens, Treat, 
                             Trial.ID, Weighted=TRUE, Alpha=.05){
  
  Surr <- Dataset[,paste(substitute(Surr))]
  SurrCens <- Dataset[,paste(substitute(SurrCens))] 
  True <- Dataset[,paste(substitute(True))]
  TrueCens <- Dataset[,paste(substitute(TrueCens))]
  Treat <- Dataset[,paste(substitute(Treat))]
  Trial.ID <- Dataset[,paste(substitute(Trial.ID))]
  All_data <- data.frame(cbind(Surr, SurrCens, True, TrueCens, Treat, Trial.ID))
  names(All_data) <- c("Surr", "SurrCens", "True", "TrueCens", "Treat", "Trial.ID")
  
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
  colnames(Trial.R2) <- c("R2.ht", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(Trial.R2) <- c(" ") 
  
  # Trial R
  Trial.R.value <- sqrt(as.numeric(summary(Results.Stage.2)[c("r.squared")]))
  Z <- .5*log((1+Trial.R.value)/(1-Trial.R.value)) 
  Trial.R.lb <- max(0, (exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))-1)/(exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))+1))
  Trial.R.ub <- min(1, (exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))-1)/(exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))+1))
  Trial.R.sd <- sqrt((1-Trial.R.value**2)/(N.trial-2))
  Trial.R <- data.frame(cbind(Trial.R.value, Trial.R.sd , Trial.R.lb, Trial.R.ub))
  colnames(Trial.R) <- c("R.ht", "Standard Error", "CI lower limit", "CI upper limit")
  row.names(Trial.R) <- c(" ") 
  

  fit <- 
    list(Data.Analyze=Data.Analyze, Results.Stage.1=Results.Stage.1, 
         Results.Stage.2=Results.Stage.2, R2.ht=Trial.R2, R.ht=Trial.R, Call=match.call())
  
  class(fit) <- "TwoStageSurvSurv"
  fit
  
}


