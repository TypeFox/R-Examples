BifixedContCont <- function(Dataset, Surr, True, Treat, Trial.ID, Pat.ID, Model=c("Full"), 
                            Weighted=TRUE, Min.Trial.Size=2, Alpha=.05){          
  
  if ((Model==c("Full") | Model==c("Reduced") | Model==c("SemiReduced"))==FALSE) {stop ("The specification of the Model=c(\"...\") argument of the call is incorrect. Use either Model=c(\"Full\"), Model=c(\"Reduced\"), or Model=c(\"SemiReduced\").")}     
  
  Surr <- Dataset[,paste(substitute(Surr))]
  True <- Dataset[,paste(substitute(True))]
  Treat <- Dataset[,paste(substitute(Treat))]
  Trial.ID <- Dataset[,paste(substitute(Trial.ID))]
  Pat.ID <- Dataset[,paste(substitute(Pat.ID))]
  
  # call Data.Processing
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
    outcomes.S.T <- cbind(dataS$outcome, dataT$outcome)
    colnames(outcomes.S.T) <- c("Surrogate", "True")
    Model.stage.1 <- lm(outcomes.S.T ~ -1 + as.factor(dataS$Trial.ID) + as.factor(dataS$Trial.ID):dataS$Treat, data=dataS)
    Intercept.S <- Model.stage.1$coefficients[1:N.trial]
    Intercept.T <- Model.stage.1$coefficients[,2][1:N.trial]
    Treatment.S <- Model.stage.1$coefficients[(N.trial+1):(2*N.trial)]
    Treatment.T <- Model.stage.1$coefficients[,2][(N.trial+1):(2*N.trial)]
    Residuals.Model.S <- residuals(Model.stage.1)[,1]
    Residuals.Model.T <- residuals(Model.stage.1)[,2]
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
    outcomes.S.T <- cbind(dataS$outcome, dataT$outcome)
    colnames(outcomes.S.T) <- c("Surrogate", "True")
    Model.stage.1 <- lm(outcomes.S.T ~ as.factor(dataS$Trial.ID):dataS$Treat, data=dataS)
    Treatment.S <- Model.stage.1$coefficients[2:(N.trial+1)]
    Treatment.T <- Model.stage.1$coefficients[,2][2:(N.trial+1)]
    Residuals.Model.S <- residuals(Model.stage.1)[,1]
    Residuals.Model.T <- residuals(Model.stage.1)[,2]
    Results.Stage.1 <- data.frame(Obs.per.trial$Trial, Obs.per.trial$Obs.per.trial, Treatment.S, Treatment.T)
    colnames(Results.Stage.1) <- c(NULL, "Trial", "Obs.per.trial", "Treatment.S", "Treatment.T")
    rownames(Results.Stage.1) <- NULL 
    D.equiv <- var(Results.Stage.1[,3:4])
    Residuals.Stage.1 <- data.frame(Residuals.Model.T, Residuals.Model.S)
    rownames(Residuals.Stage.1) <- NULL
    Residuals.Stage.1 <- cbind(wide$Pat.ID, Residuals.Stage.1)
    colnames(Residuals.Stage.1) <- c("Pat.ID", "Residuals.Model.T", "Residuals.Model.S") 
  }
  
  # Stage 2
  if (Model==c("Full")){
    if (Weighted == FALSE){
      Results.Stage.2 <- lm(Treatment.T ~ Intercept.S + Treatment.S, data=Results.Stage.1)}
    if (Weighted == TRUE){
      Results.Stage.2 <- lm(Treatment.T ~ Intercept.S + Treatment.S, data=Results.Stage.1, weights=Results.Stage.1$Obs.per.trial)}
  }
  if (Model==c("Reduced") | Model==c("SemiReduced")){
    if (Weighted == FALSE){
      Results.Stage.2 <- lm(Treatment.T ~ Treatment.S, data=Results.Stage.1)}
    if (Weighted == TRUE){
      Results.Stage.2 <- lm(Treatment.T ~ Treatment.S, data=Results.Stage.1, weights=Results.Stage.1$Obs.per.trial)}
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
  
  # Individual-Level Surrogacy and CI
  if (Model==c("Full") | Model==c("SemiReduced")) {
    Model.For.R2indiv <- nlme::gls(outcome ~ -1 + as.factor(Trial.ID):Treat:as.factor(endpoint) + as.factor(Trial.ID):as.factor(endpoint),
                             correlation = nlme::corSymm(form=~1|as.factor(Trial.ID)/as.factor(Pat.ID)), 
                             weights=nlme::varIdent(form=~1|endpoint), data=Data.analyze) }
  if (Model==c("Reduced")){
    Model.For.R2indiv <- nlme::gls(outcome ~ -1 + as.factor(Trial.ID):Treat:as.factor(endpoint) + as.factor(endpoint), 
                             correlation = nlme::corSymm(form=~1|as.factor(Trial.ID)/as.factor(Pat.ID)), 
                             weights=nlme::varIdent(form=~1|endpoint), data=Data.analyze) }
  
  cors <- nlme::corMatrix(Model.For.R2indiv$modelStruct$corStruct)[[1]]
  varStruct <- capture.output(Model.For.R2indiv$modelStruct$varStruct)[3]
  varStruct <- cbind (as.numeric(unique(strsplit(varStruct, " ")[[1]])[1]), as.numeric(unique(strsplit(varStruct, " ")[[1]])[2]))
  vars <- as.numeric((varStruct**2) * (summary(Model.For.R2indiv)$sigma)**2)
  covs <- outer(vars, vars, function(x,y) sqrt(x)*sqrt(y))
  VarCovarResid <- cors * covs
  rownames(VarCovarResid) <- colnames(VarCovarResid) <- c("Surr", "True")
  R2ind <- (VarCovarResid[2,1]**2)/(VarCovarResid[1,1]*VarCovarResid[2,2])
  R2ind.sd <- sqrt((4*R2ind*((1-R2ind)**2))/(N.total-3))
  R2ind.lb <- max(0,R2ind + qnorm(Alpha/2)*R2ind.sd)
  R2ind.ub <- min(1,R2ind + qnorm(1-Alpha/2)*R2ind.sd)
  Indiv.R2 <- data.frame(cbind(R2ind, R2ind.sd, R2ind.lb, R2ind.ub))
  colnames(Indiv.R2) <- c("R2 Indiv", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(Indiv.R2) <- c(" ")
  Min.Eigen.VarCovarResid <- min(eigen(VarCovarResid)$values)    # lowest eigenvalue
  if (Min.Eigen.VarCovarResid <= 0) warning(paste("The R-square Individual estimate may be invalid, because its calculation is based on a non-positive definite covariance matrix. "))
  
  # Rind
  Rind <- sqrt((VarCovarResid[2,1]**2)/(VarCovarResid[1,1]*VarCovarResid[2,2]))
  Z <- .5*log((1+Rind)/(1-Rind))
  Indiv.R.lb <- max(0, (exp(2*(Z-(qnorm(1-.05/2)*sqrt(1/(N.trial-3)))))-1)/(exp(2*(Z-(qnorm(1-.05/2)*sqrt(1/(N.trial-3)))))+1))
  Indiv.R.ub <- min(1, (exp(2*(Z+(qnorm(1-.05/2)*sqrt(1/(N.trial-3)))))-1)/(exp(2*(Z+(qnorm(1-.05/2)*sqrt(1/(N.trial-3)))))+1))
  Indiv.R.sd <- sqrt((1-Rind**2)/(N.total-2))
  Indiv.R <- data.frame(cbind(Rind, Indiv.R.sd, Indiv.R.lb, Indiv.R.ub))
  colnames(Indiv.R) <- c("R Indiv", "Standard Error", "CI lower limit", "CI upper limit")
  row.names(Indiv.R) <- c(" ") 
  
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
         D.Equiv=D.equiv, Sigma=VarCovarResid, Call=match.call())   
  
  class(fit) <- "BifixedContCont"
  fit
  
}
