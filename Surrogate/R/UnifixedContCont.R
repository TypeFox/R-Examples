UnifixedContCont <- function(Dataset, Surr, True, Treat, Trial.ID, Pat.ID, Model=c("Full"), Weighted=TRUE, Min.Trial.Size=2, 
                     Alpha=.05, Number.Bootstraps=500, Seed=sample(1:1000, size=1)){
  
  if ((Model==c("Full") | Model==c("Reduced") | Model==c("SemiReduced"))==FALSE) {stop ("The specification of the Model=c(\"...\") argument of the call is incorrect. Use either Model=c(\"Full\"), Model=c(\"Reduced\"), or Model=c(\"SemiReduced\").")}     
  Surr <- Dataset[,paste(substitute(Surr))]
  True <- Dataset[,paste(substitute(True))]
  Treat <- Dataset[,paste(substitute(Treat))]
  Trial.ID <- Dataset[,paste(substitute(Trial.ID))]
  Pat.ID <- Dataset[,paste(substitute(Pat.ID))]
  
  # Call Data.Processing
  Data.Proc <- .Data.Processing(Dataset=Dataset, Surr=Surr, True=True, Treat=Treat, Trial.ID=Trial.ID, Pat.ID=Pat.ID, Min.Trial.Size=Min.Trial.Size)
  wide <- Data.Proc$wide
  dataS <- Data.Proc$dataS
  dataT <- Data.Proc$dataT
  Data.analyze <- Data.Proc$Data.analyze
  N.total <- Data.Proc$N.total
  N.trial <- Data.Proc$N.trial
  Obs.per.trial <- Data.Proc$Obs.per.trial
  
  # stage 1
  if (Model==c("Full")|Model==c("SemiReduced")){
    Model.S <- lm(outcome ~ -1 + as.factor(Trial.ID) + as.factor(Trial.ID):Treat, data=dataS)      
    Model.T <- lm(outcome ~ -1 + as.factor(Trial.ID) + as.factor(Trial.ID):Treat, data=dataT)
    Intercept.S <- coef(Model.S)[1:N.trial] 
    Treatment.S <- coef(Model.S)[(1+N.trial):(2*N.trial)]
    Intercept.T <- coef(Model.T)[1:N.trial]
    Treatment.T <- coef(Model.T)[(1+N.trial):(2*N.trial)]
    Residuals.Model.S <- Model.S$residuals
    Residuals.Model.T <- Model.T$residuals
    Results.Stage.1 <- data.frame(Obs.per.trial$Trial, Obs.per.trial$Obs.per.trial, Intercept.S, Intercept.T, Treatment.S, Treatment.T)         
    colnames(Results.Stage.1) <- c(NULL, "Trial", "Obs.per.trial", "Intercept.S", "Intercept.T", "Treatment.S", "Treatment.T")
    rownames(Results.Stage.1) <- NULL 
    D.equiv <- var(Results.Stage.1[,3:6])
    Residuals.Stage.1 <- data.frame(Residuals.Model.T, Residuals.Model.S) 
    rownames(Residuals.Stage.1) <- NULL
    Residuals.Stage.1 <- cbind(wide$Pat.ID, Residuals.Stage.1)
    colnames(Residuals.Stage.1) <- c("Pat.ID", "Residuals.Model.T", "Residuals.Model.S")  
  }
  
  if (Model==c("Reduced")){  
    Model.S <- lm(outcome ~ as.factor(Trial.ID):Treat, data=dataS)      
    Model.T <- lm(outcome ~ as.factor(Trial.ID):Treat, data=dataT)
    Treatment.S <- coef(Model.S)[2:(N.trial+1)]
    Treatment.T <- coef(Model.T)[2:(N.trial+1)]
    Residuals.Model.S <- Model.S$residuals
    Residuals.Model.T <- Model.T$residuals
    Results.Stage.1 <- data.frame(Obs.per.trial$Trial, Obs.per.trial$Obs.per.trial, Treatment.S, Treatment.T)
    colnames(Results.Stage.1) <- c(NULL, "Trial", "Obs.per.trial", "Treatment.S", "Treatment.T")
    rownames(Results.Stage.1) <- NULL
    D.equiv <- var(Results.Stage.1[,3:4])     # "Treatment.S", "Treatment.T
    Residuals.Stage.1 <- data.frame(Residuals.Model.T, Residuals.Model.S)
    rownames(Residuals.Stage.1) <- NULL
    Residuals.Stage.1 <- cbind(wide$Pat.ID, Residuals.Stage.1)
    colnames(Residuals.Stage.1) <- c("Pat.ID", "Residuals.Model.T", "Residuals.Model.S") 
  }
  
  # stage 2
  if (Weighted==FALSE){
    if (Model==c("Full")){
      Results.Stage.2 <- lm(Results.Stage.1$Treatment.T ~ Results.Stage.1$Intercept.S + Results.Stage.1$Treatment.S)
    }
    if (Model==c("Reduced") | Model==c("SemiReduced")){
      Results.Stage.2 <- lm(Results.Stage.1$Treatment.T ~ Results.Stage.1$Treatment.S)
    }
  }
  if (Weighted==TRUE){
    if (Model==c("Full")){
      Results.Stage.2 <- lm(Results.Stage.1$Treatment.T ~ Results.Stage.1$Intercept.S + Results.Stage.1$Treatment.S, weights=Results.Stage.1$Obs.per.trial) 
    }
    if (Model==c("Reduced") | Model==c("SemiReduced")){
      Results.Stage.2 <- lm(Results.Stage.1$Treatment.T ~ Results.Stage.1$Treatment.S, weights=Results.Stage.1$Obs.per.trial) 
    }
  }
  
  # Trial R2
  Trial.R2.value <- as.numeric(summary(Results.Stage.2)[c("r.squared")])
  Trial.R2.sd <- sqrt((4*Trial.R2.value*(1-Trial.R2.value)^2)/(N.trial-3))
  Trial.R2.lb <- max(0, Trial.R2.value + qnorm(Alpha/2) *(Trial.R2.sd))
  Trial.R2.ub <- min(1, Trial.R2.value + qnorm(1-Alpha/2)*(Trial.R2.sd))
  Trial.R2 <- data.frame(cbind(Trial.R2.value, Trial.R2.sd, Trial.R2.lb, Trial.R2.ub))
  colnames(Trial.R2) <- c("R2 Trial", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(Trial.R2) <- c(" ") 
  
  # Trial R
  Trial.R.value <- sqrt(as.numeric(summary(Results.Stage.2)[c("r.squared")]))
  Z <- .5*log((1+Trial.R.value)/(1-Trial.R.value)) 
  Trial.R.lb <- max(0, (exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))-1)/(exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))+1))
  Trial.R.ub <- min(1, (exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))-1)/(exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))+1))
  Trial.R.sd <- sqrt((1-Trial.R.value**2)/(N.trial-2))
  Trial.R <- data.frame(cbind(Trial.R.value, Trial.R.sd , Trial.R.lb, Trial.R.ub))
  colnames(Trial.R) <- c("R Trial", "Standard Error", "CI lower limit", "CI upper limit")
  row.names(Trial.R) <- c(" ") 
  
  # R2 indiv
  R2ind <- (cor(Residuals.Model.T, Residuals.Model.S))**2
  Boot.r <-rep(0, Number.Bootstraps)
  for (j in 1:Number.Bootstraps){
    obs <- c(1:N.total)
    set.seed(Seed)
    Indicator <- sample(obs, N.total, replace=TRUE)  
    Seed <- Seed + 1
    Sample.boot.S <- data.frame(dataS[Indicator,])
    Sample.boot.T <- data.frame(dataT[Indicator,])
    Sample.boot.S <- na.exclude(Sample.boot.S[order(Sample.boot.S$Pat.ID),])   
    Sample.boot.T <- na.exclude(Sample.boot.T[order(Sample.boot.T$Pat.ID),])
    if (Model==c("Full") | Model==c("SemiReduced")){
      Boot.model.S <- lm(outcome ~ -1 + as.factor(Trial.ID) + as.factor(Trial.ID):Treat, data=Sample.boot.S)
      Boot.model.T <- lm(outcome ~ -1 + as.factor(Trial.ID) + as.factor(Trial.ID):Treat, data=Sample.boot.T)
    }
    if (Model==c("Reduced") ){
      Boot.model.S <- lm(outcome ~ as.factor(Trial.ID):Treat, data=Sample.boot.S)
      Boot.model.T <- lm(outcome ~ as.factor(Trial.ID):Treat, data=Sample.boot.T)
    }
    
    Res.Boot.model.S<-residuals(Boot.model.S, type='response')
    Res.Boot.model.T<-residuals(Boot.model.T, type='response')
    Boot.r[j] <- (cor(Res.Boot.model.S,Res.Boot.model.T))^2
  }
  
  # nonparametric bootstrap normal confidence interval
  Var.Boot.r <- var(Boot.r)
  Indiv.R2.lb <- max(0, R2ind + qnorm(Alpha/2)*sqrt(Var.Boot.r))
  Indiv.R2.ub <- min(1, R2ind - qnorm(Alpha/2)*sqrt(Var.Boot.r))
  Indiv.R2 <- data.frame(cbind(R2ind, sqrt(Var.Boot.r), Indiv.R2.lb, Indiv.R2.ub))
  colnames(Indiv.R2) <- c("R2 Indiv", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(Indiv.R2) <- c(" ")
  
  # Rind
  Rind <- cor(Residuals.Model.T, Residuals.Model.S)
  Boot.r <-rep(0, Number.Bootstraps)
  for (j in 1:Number.Bootstraps){
    obs <- c(1:N.total)
    Indicator <- sample(obs, N.total, replace=TRUE) 
    Sample.boot.S <- data.frame(dataS[Indicator,])
    Sample.boot.T <- data.frame(dataT[Indicator,])
    Sample.boot.S <- na.exclude(Sample.boot.S[order(Sample.boot.S$Pat.ID),])      
    Sample.boot.T <- na.exclude(Sample.boot.T[order(Sample.boot.T$Pat.ID),])
    Boot.model.S <- lm(outcome ~ as.factor(Trial.ID):Treat, data=Sample.boot.S)
    Boot.model.T <- lm(outcome ~ as.factor(Trial.ID):Treat, data=Sample.boot.T)
    Res.Boot.model.S<-residuals(Boot.model.S, type='response')
    Res.Boot.model.T<-residuals(Boot.model.T, type='response')
    Boot.r[j] <- (cor(Res.Boot.model.S,Res.Boot.model.T))^2
  }
  
  # nonparametric bootstrap normal confidence interval
  Var.Boot.r <- var(Boot.r)
  Indiv.R.lb <- max(0, Rind + qnorm(Alpha/2)*sqrt(Var.Boot.r))
  Indiv.R.ub <- min(1, Rind - qnorm(Alpha/2)*sqrt(Var.Boot.r))
  Indiv.R <- data.frame(cbind(Rind, sqrt(Var.Boot.r), Indiv.R.lb, Indiv.R.ub))
  colnames(Indiv.R) <- c("R Indiv", "Standard Error", "CI lower limit", "CI upper limit") 
  rownames(Indiv.R) <- c(" ") 
  
  NoTreat <- wide[wide$Treat!=1,]
  Treat <- wide[wide$Treat==1,]
  T0S0 <- cor(NoTreat$Surr, NoTreat$True)
  T1S1 <- cor(Treat$Surr, Treat$True)  
  Z_T0S0 <- .5*log((1+T0S0)/(1-T0S0)) 
  rho_lb <- max(0, (exp(2*(Z_T0S0-(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))-1)/(exp(2*(Z_T0S0-(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))+1))
  rho_ub <- min(1, (exp(2*(Z_T0S0+(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))-1)/(exp(2*(Z_T0S0+(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))+1))
  rho_sd <- sqrt((1-T0S0**2)/(N.total-2))
  rho_results_T0S0 <- data.frame(cbind(T0S0, rho_sd , rho_lb, rho_ub))
  colnames(rho_results_T0S0) <- c("Estimate", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(rho_results_T0S0) <- c(" ")
  Z_T1S1 <- .5*log((1+T1S1)/(1-T1S1)) 
  rho_lb <- max(0, (exp(2*(Z_T1S1-(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))-1)/(exp(2*(Z_T1S1-(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))+1))
  rho_ub <- min(1, (exp(2*(Z_T1S1+(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))-1)/(exp(2*(Z_T1S1+(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))+1))
  rho_sd <- sqrt((1-T1S1**2)/(N.total-2))
  rho_results_T1S1 <- data.frame(cbind(T1S1, rho_sd , rho_lb, rho_ub))
  colnames(rho_results_T1S1) <- c("Estimate", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(rho_results_T1S1) <- c(" ")
  Cor.Endpoints <- data.frame(rbind(rho_results_T0S0, rho_results_T1S1))
  rownames(Cor.Endpoints) <- c("r_T0S0", "r_T1S1")
  colnames(Cor.Endpoints) <- c("Estimate", "Standard Error", "CI lower limit", "CI upper limit")
  
  fit <- 
    list(Data.Analyze=wide, Obs.Per.Trial=Obs.per.trial, Results.Stage.1=Results.Stage.1, Residuals.Stage.1=Residuals.Stage.1, 
         Results.Stage.2=Results.Stage.2, Trial.R2=Trial.R2, Indiv.R2=Indiv.R2, Trial.R=Trial.R, Indiv.R=Indiv.R, Cor.Endpoints=Cor.Endpoints, 
         D.Equiv=D.equiv, Call=match.call())
  
  class(fit) <- "UnifixedContCont"
  fit
  
}
