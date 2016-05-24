BimixedContCont <- function(Dataset, Surr, True, Treat, Trial.ID, Pat.ID, Model=c("Full"), 
                    Min.Trial.Size=2, Alpha=.05, ...){
  
  if ((Model==c("Full") | Model==c("Reduced"))==FALSE) {stop ("The specification of the Model=c(\"...\") argument of the call is incorrect. Use either Model=c(\"Full\") or Model=c(\"Reduced\").")}     
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
  
  if (Model==c("Full")){
    
    Model.S.and.T <- lmer(outcome ~ -1 + as.factor(endpoint):Treat + as.factor(endpoint) + (-1+endpoint+endpoint:Treat|Trial.ID), data=Data.analyze, ...)   #
    
    Fixed.effect.pars.S <- rbind(summary(Model.S.and.T)$coefficients[1], summary(Model.S.and.T)$coefficients[3]) 
    Fixed.effect.pars.T <- rbind(summary(Model.S.and.T)$coefficients[2], summary(Model.S.and.T)$coefficients[4]) 
    Fixed.Effect.Pars <- data.frame(rbind(Fixed.effect.pars.S, Fixed.effect.pars.T))
    rownames(Fixed.Effect.Pars) <- c("Intercept.S", "Treatment.S", "Intercept.T", "Treatment.T")  
    colnames(Fixed.Effect.Pars) <- c(" ")  
    
    Random.effect.pars.S <- cbind(data.frame(ranef(Model.S.and.T)$Trial.ID)[,1], data.frame(ranef(Model.S.and.T)$Trial.ID)[,3])
    Random.effect.pars.T <- cbind(data.frame(ranef(Model.S.and.T)$Trial.ID)[,2], data.frame(ranef(Model.S.and.T)$Trial.ID)[,4])
    Random.Effect.Pars <- data.frame(cbind(Random.effect.pars.S, Random.effect.pars.T))
    colnames(Random.Effect.Pars) <- c("Intercept.S", "Treatment.S", "Intercept.T", "Treatment.T")
    
    Residuals.S.and.T <- residuals(Model.S.and.T, type='response')
    Residuals.S <- Residuals.S.and.T[seq(from=1, to=length(Residuals.S.and.T), by=2)]
    names(Residuals.S) <- NULL
    Residuals.T <- Residuals.S.and.T[seq(from=2, to=length(Residuals.S.and.T), by=2)]
    names(Residuals.T) <- NULL
    Residuals <- data.frame(cbind(wide$Pat.ID, Residuals.S, Residuals.T))
    colnames(Residuals) <- c("Pat.ID", "Residuals.S", "Residuals.T") 
    rownames(Residuals) <- NULL
    Residuals <- Residuals[order(Residuals$Pat.ID),]
    
    Intercept.S <- coef(Model.S.and.T)$Trial.ID[,1]+coef(Model.S.and.T)$Trial.ID[,5] 
    Intercept.T <- coef(Model.S.and.T)$Trial.ID[,2]+coef(Model.S.and.T)$Trial.ID[,6]  
    Treatment.S <- coef(Model.S.and.T)$Trial.ID[,3]+coef(Model.S.and.T)$Trial.ID[,7] 
    Treatment.T <- coef(Model.S.and.T)$Trial.ID[,4]+coef(Model.S.and.T)$Trial.ID[,8] 
    Trial.Spec.Results <- cbind(Obs.per.trial$Trial, Obs.per.trial$Obs.per.trial, Intercept.S, Intercept.T, Treatment.S, Treatment.T)
    colnames(Trial.Spec.Results) <- c(NULL, "Trial", "Obs.per.trial", "Intercept.S", "Intercept.T", "Treatment.S", "Treatment.T")
    rownames(Trial.Spec.Results) <- NULL 
  }
  
  if (Model==c("Reduced")){
    
    Model.S.and.T <- lmer(outcome ~ -1 + as.factor(endpoint):Treat + as.factor(endpoint) + (-1+endpoint:Treat|Trial.ID), data=Data.analyze, ...)   #
        
    Fixed.effect.pars.S <- rbind(summary(Model.S.and.T)$coefficients[1], summary(Model.S.and.T)$coefficients[3])
    Fixed.effect.pars.T <- rbind(summary(Model.S.and.T)$coefficients[2], summary(Model.S.and.T)$coefficients[4]) 
    Fixed.Effect.Pars <- rbind(Fixed.effect.pars.S, Fixed.effect.pars.T)
    rownames(Fixed.Effect.Pars) <- c("Intercept.S", "Treatment.S", "Intercept.T", "Treatment.T")
    colnames(Fixed.Effect.Pars) <- c(" ")
    
    Random.Effect.Pars <- data.frame(ranef(Model.S.and.T)$Trial.ID)
    colnames(Random.Effect.Pars) <- c("Treatment.S", "Treatment.T")
    
    Residuals.S.and.T <- residuals(Model.S.and.T, type='response')
    Residuals.S <- Residuals.S.and.T[seq(from=1, to=length(Residuals.S.and.T), by=2)]
    names(Residuals.S) <- NULL
    Residuals.T <- Residuals.S.and.T[seq(from=2, to=length(Residuals.S.and.T), by=2)]
    names(Residuals.T) <- NULL
    Residuals <- data.frame(cbind(wide$Pat.ID, Residuals.S, Residuals.T))
    colnames(Residuals) <- c("Pat.ID", "Residuals.S", "Residuals.T") 
    rownames(Residuals) <- NULL
    Residuals <- Residuals[order(Residuals$Pat.ID),]
    
    Treatment.S <- coef(Model.S.and.T)$Trial.ID[,1]+coef(Model.S.and.T)$Trial.ID[,5]
    Treatment.T <- coef(Model.S.and.T)$Trial.ID[,2]+coef(Model.S.and.T)$Trial.ID[,6]
    Trial.Spec.Results <- cbind(Obs.per.trial$Trial, Obs.per.trial$Obs.per.trial, Treatment.S, Treatment.T)
    colnames(Trial.Spec.Results) <- c(NULL, "Trial", "Obs.per.trial", "Treatment.S", "Treatment.T")
    rownames(Trial.Spec.Results) <- NULL 
  }
  
  
  # Trial-level surrogacy estimates
  if (Model==c("Full")){
    
    D <- matrix(summary(Model.S.and.T)$varcor$Trial.ID[1:16], ncol=4)
    rownames(D) <- colnames(D) <- c("Intercept.S", "Intercept.T", "Treatment.S", "Treatment.T")
    
    Min.Eigen.D <- min(eigen(D)$values)
    if (Min.Eigen.D <= 0) warning(paste("The R-square Trial estimate may be invalid, because its calculation is based on a non-positive definite covariance matrix"))
    singular  <- svd(D)$d
    Cond.Number.D.Matrix <- max(singular)/min(singular)
    if (Cond.Number.D.Matrix > 100) warning(paste("The R-square Trial estimate may not be thrustworthy, as the conditioning number of the D matrix is high and equals",
                                                  deparse(Cond.Number.D.Matrix)))  
    A <- matrix(c(D[4,1], D[4,3]), nrow=2, ncol=1)
    B <- matrix(c(D[1,1], D[3,1], D[3,1], D[3,3]), nrow=2, ncol=2)
    C <- as.matrix(D[4,4])
    Trial.R2.value <- (t(A) %*% solve(B) %*% A)/C
  }
  
  if (Model==c("Reduced")){
    
    D <- matrix(summary(Model.S.and.T)$varcor$Trial.ID[1:4], ncol=2)
    rownames(D) <- colnames(D) <- c("Treatment.S", "Treatment.T")
    Min.Eigen.D <- min(eigen(D)$values)
    if (Min.Eigen.D <= 0) warning(paste("The R-square Trial estimate may be invalid, because its calculation is based on a non-positive definite covariance matrix"))
    singular  <- svd(D)$d
    Cond.Number.D.Matrix <- max(singular)/min(singular)
    if (Cond.Number.D.Matrix > 100) warning(paste("The R-square Trial estimate may not be thrustworthy, as the conditioning number of the D matrix is high and equals",
                                                  deparse(Cond.Number.D.Matrix)))
    Trial.R2.value <- (D[2,1]**2)/(D[1,1]*D[2,2])
  }
  
  # R2 trial
  Trial.R2.sd <- sqrt((4*Trial.R2.value*(1-Trial.R2.value)^2)/(N.trial-3))
  Trial.R2.lb <- max(0, Trial.R2.value + qnorm(Alpha/2) *(Trial.R2.sd))
  Trial.R2.ub <- min(1, Trial.R2.value + qnorm(1-Alpha/2)*(Trial.R2.sd))
  Trial.R2 <- data.frame(cbind(Trial.R2.value, Trial.R2.sd, Trial.R2.lb, Trial.R2.ub))
  colnames(Trial.R2) <- c("R2 Trial", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(Trial.R2) <- c(" ")
  
  # Rtrial
  Trial.R.value <- sqrt(Trial.R2.value)
  Z <- .5*log((1+Trial.R.value)/(1-Trial.R.value)) 
  Trial.R.lb <- max(0, (exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))-1)/(exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))+1))
  Trial.R.ub <- min(1, (exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))-1)/(exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(N.trial-3)))))+1))
  Trial.R.sd <- sqrt((1-Trial.R.value**2)/(N.trial-2))
  Trial.R <- data.frame(cbind(Trial.R.value, Trial.R.sd, Trial.R.lb, Trial.R.ub))
  colnames(Trial.R) <- c("R Trial", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(Trial.R) <- c(" ")
  
  # R2 indiv
  
  if (Model==c("Full")){
    
    rm(Model.S.and.T)
    Model.S.and.T <- lme(outcome~ -1 + as.factor(endpoint):Treat + as.factor(endpoint), 
                     random=~ -1 + as.factor(endpoint) + as.factor(endpoint):Treat|as.factor(Trial.ID),
                     correlation = nlme::corSymm(form=~1|as.factor(Trial.ID)/as.factor(Pat.ID)), data=Data.analyze,
                     weights=nlme::varIdent(form=~1|endpoint), control = list(msVerbose = FALSE, optimizer = "nlm", niterEM = 25, msMaxIter=500))  
   }
  
  if (Model==c("Reduced")){
    rm(Model.S.and.T)
    Model.S.and.T <- lme(outcome~ -1 + as.factor(endpoint):Treat + as.factor(endpoint), 
                         random=~ -1 + as.factor(endpoint):Treat|as.factor(Trial.ID),
                         correlation = nlme::corSymm(form=~1|as.factor(Trial.ID)/as.factor(Pat.ID)), data=Data.analyze,
                         weights=nlme::varIdent(form=~1|endpoint), control = list(msVerbose = FALSE, optimizer = "nlm", niterEM = 25, msMaxIter=500))  
   }
  
  cors <- nlme::corMatrix(Model.S.and.T$modelStruct$corStruct)[[1]]
  varStruct <- capture.output(Model.S.and.T$modelStruct$varStruct)[3]
  varStruct <- cbind(as.numeric(unique(strsplit(varStruct, " ")[[1]])[1]), as.numeric(unique(strsplit(varStruct, " ")[[1]])[2]))
  vars <- as.numeric((varStruct**2) * (summary(Model.S.and.T)$sigma)**2)
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
  
  # R ind
  Rind <- sqrt((VarCovarResid[2,1]**2)/(VarCovarResid[1,1]*VarCovarResid[2,2]))
  Z <- .5*log((1+Rind)/(1-Rind))
  Indiv.R.lb <- max(0, (exp(2*(Z-(qnorm(1-.05/2)*sqrt(1/(N.trial-3)))))-1)/(exp(2*(Z-(qnorm(1-.05/2)*sqrt(1/(N.trial-3)))))+1))
  Indiv.R.ub <- min(1, (exp(2*(Z+(qnorm(1-.05/2)*sqrt(1/(N.trial-3)))))-1)/(exp(2*(Z+(qnorm(1-.05/2)*sqrt(1/(N.trial-3)))))+1))
  Indiv.R.sd <- sqrt((1-Rind**2)/(N.total-2))
  Indiv.R <- data.frame(cbind(Rind, Indiv.R.sd, Indiv.R.lb, Indiv.R.ub))
  colnames(Indiv.R) <- c("R Indiv", "Standard Error", "CI lower limit", "CI upper limit")
  row.names(Indiv.R) <- c(" ") 
  
  Min.Eigen.VarCovarResid <- min(eigen(VarCovarResid)$values)    # lowest eigenvalue
  if (Min.Eigen.VarCovarResid <= 0) warning(paste("The R-square Individual estimate may be invalid, because its calculation is based on a non-positive definite covariance matrix"))
  Singular <- svd(VarCovarResid)$d
  Cond.Number.VarCovarResid <- max(Singular)/min(Singular)
  if (Cond.Number.VarCovarResid > 100) warning(paste("The R-square Individual estimate may not be thrustworthy, as the conditioning number of the Sigma matrix is high and equals",
                                                     deparse(Cond.Number.VarCovarResid)))
  
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
    list(Data.Analyze=wide, Obs.Per.Trial=Obs.per.trial, Trial.Spec.Results=Trial.Spec.Results, Residuals=Residuals, 
         Fixed.Effect.Pars=Fixed.Effect.Pars, Random.Effect.Pars=Random.Effect.Pars, Trial.R2=Trial.R2, Trial.R=Trial.R, 
         Indiv.R2=Indiv.R2, Indiv.R=Indiv.R, Cor.Endpoints=Cor.Endpoints, D=D, Sigma=VarCovarResid, Call=match.call())
  
  class(fit) <- "BimixedContCont"
  
  fit
  
}