UnimixedContCont <- function(Dataset, Surr, True, Treat, Trial.ID, Pat.ID, Model=c("Full"), 
                     Weighted=TRUE, Min.Trial.Size=2, Alpha=.05, Number.Bootstraps=500, 
                     Seed=sample(1:1000, size=1), ...){
  
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
  
  Control=list(msMaxIter=500)
  
  # Stage 1
  # fit univariate mixed-effect models for S and T
  if (Model==c("Full")|Model==c("SemiReduced")){
    
    Model.S <- lmer(outcome ~ Treat+(1+Treat|Trial.ID), data=dataS, ...)  
    Model.T <- lmer(outcome ~ Treat+(1+Treat|Trial.ID), data=dataT, ...)
    
    Intercept.S <- coef(Model.S)$Trial.ID[,1] #coefficients
    Treatment.S <- coef(Model.S)$Trial.ID[,2]
    Intercept.T <- coef(Model.T)$Trial.ID[,1]
    Treatment.T <- coef(Model.T)$Trial.ID[,2]
    Results.Stage.1 <- data.frame(Obs.per.trial$Trial, Obs.per.trial$Obs.per.trial, Intercept.S, Intercept.T, Treatment.S, Treatment.T)         
    colnames(Results.Stage.1) <- c(NULL, "Trial", "Obs.per.trial", "Intercept.S", "Intercept.T", "Treatment.S", "Treatment.T")
    rownames(Results.Stage.1) <- NULL
    D.equiv <- var(Results.Stage.1[,3:6])
    Residuals.Model.S <- residuals(Model.S, type='response')
    Residuals.Model.T <- residuals(Model.T, type='response')
    Residuals.Stage.1 <- cbind(wide$Pat.ID, data.frame(Residuals.Model.S, Residuals.Model.T))
    colnames(Residuals.Stage.1) <- c("Pat.ID", "Residuals.Model.S", "Residuals.Model.T") 
    rownames(Residuals.Stage.1) <- NULL
    
    Fixed.effect.pars.S <- matrix(summary(Model.S)$coef[1:2], nrow=2)  # intercept S en treat S
    Fixed.effect.pars.T <- matrix(summary(Model.T)$coef[1:2], nrow=2) 
    rownames(Fixed.effect.pars.S) <- c("Intercept.S" , "Treatment.S")  
    rownames(Fixed.effect.pars.T) <- c("Intercept.T" , "Treatment.T")
    Fixed.Effect.Pars <- data.frame(rbind(Fixed.effect.pars.S, Fixed.effect.pars.T))
    colnames(Fixed.Effect.Pars) <- c(" ")
    
    Random.effect.pars.S <- data.frame(ranef(Model.S)$Trial.ID)
    Random.effect.pars.T <- data.frame(ranef(Model.T)$Trial.ID)
    colnames(Random.effect.pars.S) <- c("Intercept.S", "Treatment.S")
    colnames(Random.effect.pars.T) <- c("Intercept.S", "Treatment.S")
    Random.Effect.Pars <- cbind(Random.effect.pars.S, Random.effect.pars.T)
  }
  
    
  if (Model==c("Reduced")){
    
    Model.S <- lmer(outcome ~ Treat+(-1+Treat|Trial.ID), data=dataS, ...)  
    Model.T <- lmer(outcome ~ Treat+(-1+Treat|Trial.ID), data=dataT, ...)
    
    Treatment.S <- coef(Model.S)$Trial.ID[,2]
    names(Treatment.S)<-"Treatment.S"
    Treatment.T <- coef(Model.T)$Trial.ID[,2]
    names(Treatment.T)<-"Treatment.T"
    Results.Stage.1 <- data.frame(Obs.per.trial$Trial, Obs.per.trial$Obs.per.trial, Treatment.S, Treatment.T)         
    colnames(Results.Stage.1) <- c(NULL, "Trial", "Obs.per.trial", "Treatment.S", "Treatment.T")
    rownames(Results.Stage.1) <- NULL 
    D.equiv <- var(Results.Stage.1[,3:4])
    
    Residuals.Model.S <- residuals(Model.S, type='response')
    Residuals.Model.T <- residuals(Model.T, type='response')
    Residuals.Stage.1 <- cbind(wide$Pat.ID, data.frame(Surr=Residuals.Model.S, True=Residuals.Model.T))
    colnames(Residuals.Stage.1) <- c("Pat.ID", "Residuals.Model.S", "Residuals.Model.T") 
    rownames(Residuals.Stage.1) <- NULL
    
    Fixed.effect.pars.S <- matrix(summary(Model.S)$coef[1:2], nrow=2) 
    rownames(Fixed.effect.pars.S)[1:2]<-c("Intercept.S", "Treatment.S")
    Fixed.effect.pars.T <- matrix(summary(Model.T)$coef[1:2], nrow=2) 
    rownames(Fixed.effect.pars.T)[1:2]<-c("Intercept.T", "Treatment.T")
    Fixed.Effect.Pars <- data.frame(rbind(Fixed.effect.pars.S, Fixed.effect.pars.T))
    colnames(Fixed.Effect.Pars) <- c(" ")
    
    Random.effect.pars.S <- data.frame(ranef(Model.S)$Trial.ID)
    colnames(Random.effect.pars.S) <- c("Treatment.S")
    Random.effect.pars.T <- data.frame(ranef(Model.T)$Trial.ID)
    colnames(Random.effect.pars.T) <- c("Treatment.T")
    Random.Effect.Pars <- cbind(Random.effect.pars.S, Random.effect.pars.T)
  }
  
  # Trial Level Surrogacy   
  
  if (Model==c("Full")){
    if (Weighted==FALSE) {Results.Stage.2 <- lm(Results.Stage.1$Treatment.T ~ Results.Stage.1$Intercept.S + Results.Stage.1$Treatment.S)}
    if (Weighted==TRUE) {Results.Stage.2 <- lm(Results.Stage.1$Treatment.T ~ Results.Stage.1$Intercept.S + Results.Stage.1$Treatment.S, weights=Results.Stage.1$Obs.per.trial)}
  }
  if (Model==c("Reduced") | Model==c("SemiReduced")){
    if (Weighted==FALSE) {Results.Stage.2 <- lm(Results.Stage.1$Treatment.T ~ Results.Stage.1$Treatment.S)}
    if (Weighted==TRUE) {Results.Stage.2 <- lm(Results.Stage.1$Treatment.T ~ Results.Stage.1$Treatment.S, weights=Results.Stage.1$Obs.per.trial)}
  }
  
  
  # R2 trial
  Trial.R2.value <- as.numeric(summary(Results.Stage.2)[c("r.squared")])
  Trial.R2.sd <- sqrt((4*Trial.R2.value*(1-Trial.R2.value)^2)/(N.trial-3))
  Trial.R2.lb <- max(0, Trial.R2.value + qnorm(Alpha/2) *(Trial.R2.sd))
  Trial.R2.ub <- min(1, Trial.R2.value + qnorm(1-Alpha/2)*(Trial.R2.sd))
  Trial.R2 <- data.frame(cbind(Trial.R2.value, Trial.R2.sd, Trial.R2.lb, Trial.R2.ub))
  colnames(Trial.R2) <- c("R2 Trial", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(Trial.R2) <- c(" ")
  
  # Rtrial
  Trial.R.value <- sqrt(as.numeric(summary(Results.Stage.2)[c("r.squared")]))
  Z <- .5*log((1+Trial.R.value)/(1-Trial.R.value)) 
  Trial.R.lb <- max(0, (exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))-1)/(exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))+1))
  Trial.R.ub <- min(1, (exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))-1)/(exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))+1))
  Trial.R.sd <- sqrt((1-Trial.R.value**2)/(N.trial-2))
  Trial.R <- data.frame(cbind(Trial.R.value, Trial.R.sd, Trial.R.lb, Trial.R.ub))
  colnames(Trial.R) <- c("R Trial", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(Trial.R) <- c(" ")
  
  # Individual Level Surrogacy
  options(warn = -1)
  Boot.r <- rep(0, Number.Bootstraps) 
  for (j in 1:Number.Bootstraps){
    obs <- c(1:N.total)
    set.seed(Seed)
    Indicator <- sample(obs, N.total, replace=TRUE)
    Seed <- Seed + 1
    Sample.boot.S <- data.frame(dataS[Indicator,])
    Sample.boot.T <- data.frame(dataT[Indicator,])
    
    if (Model==c("Full") | Model==c("SemiReduced")){                              
      Boot.model.S <- try(lmer(outcome ~ Treat+(1+Treat|Trial.ID), data=Sample.boot.S, ...), silent = FALSE)
      Boot.model.T <- try(lmer(outcome ~ Treat+(1+Treat|Trial.ID), data=Sample.boot.T, ...), silent = FALSE)
    }
    if (Model==c("Reduced")){
      Boot.model.S <- try(lmer(outcome ~ Treat+(-1+Treat|Trial.ID), data=Sample.boot.S, ...), silent = FALSE)
      Boot.model.T <- try(lmer(outcome ~ Treat+(-1+Treat|Trial.ID), data=Sample.boot.T, ...), silent = FALSE)
    }
    Res.Boot.model.S <- residuals(Boot.model.S, type='response')
    Res.Boot.model.T <- residuals(Boot.model.T, type='response')
    Boot.r[j] <- (cor(Res.Boot.model.S,Res.Boot.model.T))
  }
  Boot.r2 <- Boot.r**2
  options(warn=0)
  
  # R2 ind
  R2ind <- (cor(Residuals.Model.T, Residuals.Model.S))**2
  Var.Boot.r2 <- var(Boot.r2)
  Indiv.R2.lb <- max(0, R2ind + qnorm(Alpha/2)*sqrt(Var.Boot.r2)) 
  Indiv.R2.ub <- R2ind - qnorm(Alpha/2)*sqrt(Var.Boot.r2)
  Indiv.R2 <- data.frame(cbind(R2ind, sqrt(Var.Boot.r2), Indiv.R2.lb, Indiv.R2.ub))
  colnames(Indiv.R2) <- c("R2 Indiv", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(Indiv.R2) <- c(" ")
  
  # R ind
  Rind <- (cor(Residuals.Model.T, Residuals.Model.S))
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
         Fixed.Effect.Pars=Fixed.Effect.Pars, Random.Effect.Pars=Random.Effect.Pars, Results.Stage.2=Results.Stage.2, Trial.R2=Trial.R2, Indiv.R2=Indiv.R2, Trial.R=Trial.R, Indiv.R=Indiv.R, Cor.Endpoints=Cor.Endpoints, 
         D.Equiv=D.equiv, Call=match.call())
  
  class(fit) <- "UnimixedContCont"
  
  fit
  
}