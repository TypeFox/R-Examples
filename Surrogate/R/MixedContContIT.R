MixedContContIT <- function(Dataset, Surr, True, Treat, Trial.ID, Pat.ID,
                            Model=c("Full"), Weighted=TRUE, Min.Trial.Size=2, Alpha=.05, ...){  
  
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
  
  Trialname <- data.frame(table(Data.analyze$Trial.ID))[,1]  
  
  # R2ht and CI      
  command.Model.S <- c("Surr ~ Treat")
  command.Model.T <- c("True ~ Treat")
  
  if (Model==c("Full")|Model==c("SemiReduced")){
    command.Model.S <- paste(command.Model.S, " + (1+Treat|Trial.ID)", sep="") 
    command.Model.T <- paste(command.Model.T, " + (1+Treat|Trial.ID)", sep="")
  }
  if (Model==c("Reduced")){
    command.Model.S <- paste(command.Model.S, " + (-1+Treat|Trial.ID)", sep="") 
    command.Model.T <- paste(command.Model.T, " + (-1+Treat|Trial.ID)", sep="")
  }
  
  Model.S <- lmer(eval(parse(text=command.Model.S)), data=wide, ...) 
  Model.T <- lmer(eval(parse(text=command.Model.T)), data=wide, ...) 
  
  if (Model==c("Full")|Model==c("SemiReduced")){
    Intercept.S <- fixef(Model.S)["(Intercept)"]-ranef(Model.S)$Trial.ID["(Intercept)"][,1]  
    Intercept.T <- fixef(Model.T)["(Intercept)"]-ranef(Model.T)$Trial.ID["(Intercept)"][,1]  
  }
  
  Treatment.S <- fixef(Model.S)["Treat"]-ranef(Model.S)$Trial.ID["Treat"][,1]   
  Treatment.T <- fixef(Model.T)["Treat"]-ranef(Model.T)$Trial.ID["Treat"][,1]
  
  if (Model==c("Full")| Model==c("SemiReduced")){  
    Trial.Spec.Results <- data.frame(Obs.per.trial$Trial, Obs.per.trial$Obs.per.trial, Intercept.S, Treatment.S, Intercept.T, Treatment.T)
    colnames(Trial.Spec.Results) <- c(NULL, "Trial", "Obs.per.trial", "Intercept.S", "Treatment.S", "Intercept.T", "Treatment.T")
    rownames(Trial.Spec.Results) <- NULL 
  }
  
  if (Model==c("Reduced")){  
    Trial.Spec.Results <- data.frame(Obs.per.trial$Trial, Obs.per.trial$Obs.per.trial, Treatment.S, Treatment.T)
    colnames(Trial.Spec.Results) <- c(NULL, "Trial", "Obs.per.trial", "Treatment.S", "Treatment.T")
    rownames(Trial.Spec.Results) <- NULL
  }
  
  # residuals
  Residuals.Model.S <- residuals(Model.S) 
  Residuals.Model.T <- residuals(Model.T)
  Residuals <- data.frame(Residuals.Model.T, Residuals.Model.S) 
  rownames(Residuals) <- NULL
  Residuals <- cbind(wide$Pat.ID, Residuals)
  colnames(Residuals) <- c("Pat.ID", "Residuals.Model.T", "Residuals.Model.S")    
  
  
  if (Weighted==FALSE){
    if (Model==c("Full")){
      Model.1 <- glm(Treatment.T ~ Intercept.S + Treatment.S, family=gaussian) 
      L1 <- -2 * logLik(Model.1)[1]}
    if (Model==c("Reduced")|Model==c("SemiReduced")){
      Model.1 <- glm(Treatment.T ~ Treatment.S, family=gaussian) 
      L1 <- -2 * logLik(Model.1)[1]}
  }
  if (Weighted==TRUE){
    if (Model==c("Full")){
      Model.1 <- glm(Treatment.T ~ Intercept.S + Treatment.S, weights=Obs.per.trial$Obs.per.trial, family=gaussian) 
      L1 <- -2 * logLik(Model.1)[1]}
    if (Model==c("Reduced")|Model==c("SemiReduced")){
      Model.1 <- glm(Treatment.T ~ Treatment.S, weights=Obs.per.trial$Obs.per.trial, family=gaussian) 
      L1 <- -2 * logLik(Model.1)[1]}
  }
  
  Model.0 <- glm(Treatment.T ~ 1, family=gaussian)       
  L0 <- -2 * logLik(Model.0)[1]
  g2 <- -(L1-L0)
  R2ht.value <- 1 - exp(-g2/N.trial)
  k1 <- qchisq(Alpha, 1, g2)
  d1 <- qchisq((1-Alpha), 1, g2)
  R2ht.lb <- max(0, 1-exp(-k1/N.trial))
  R2ht.ub <- min(1, 1-exp(-d1/N.trial))
  R2ht <- data.frame(cbind(R2ht.value, R2ht.lb, R2ht.ub)) #output
  colnames(R2ht) <- c("R2ht", "CI lower limit", "CI upper limit")
  rownames(R2ht) <- c(" ")  
  
  # Individual-level
  
  # R2h.single (single-cluster based) and CI            
  command.Model.1 <- c("True ~ Treat")          
  command.Model.2 <- c("True ~ Treat + Surr")
  command.Model.1 <- paste(command.Model.1, " + (1+Treat|Trial.ID)", sep="") 
  command.Model.2 <- paste(command.Model.2, " + (1+Treat+Surr|Trial.ID)", sep="")
  
  Model.1 <- lmer(eval(parse(text=command.Model.1)), data=wide, ...)         
  Model.2 <- lmer(eval(parse(text=command.Model.2)), data=wide, ...)         
  Treatment.S <- fixef(Model.1)["Treat"]-ranef(Model.1)$Trial.ID["Treat"][,1]   
  Treatment.T <- fixef(Model.2)["Treat"]-ranef(Model.2)$Trial.ID["Treat"][,1]
  
  # R2h single-trial based
  L1 <- -2 * logLik(Model.1)[1]
  L2 <- -2 * logLik(Model.2)[1]
  g2 <- -(L2-L1)
  R2h.single.value <- 1 - exp(-g2/N.total)
  k1 <- qchisq(Alpha, 1, g2)
  d1 <- qchisq((1-Alpha), 1, g2)
  R2h.single.lb <- max(0, 1-exp(-k1/N.total)) 
  R2h.single.ub <- min(1, 1-exp(-d1/N.total))
  R2h.single <- data.frame(cbind(R2h.single.value, R2h.single.lb, R2h.single.ub))   # OUT
  colnames(R2h.single) <- c("R2h.ind", "CI lower limit", "CI upper limit")
  rownames(R2h.single) <- c(" ")
  
  # R2h.single.max and CI
  Model.0 <- glm(True ~ 1, data=wide)
  L.0 <- -2 * logLik(Model.0)[1]    
  R2h.single.max.value <- R2h.single.value/(1-(exp(-L.0/N.total)))   
  R2h.single.max.lb <- max(0, (1-exp(-k1/N.total))/(1-(exp(-L.0/N.total))))
  R2h.single.max.ub <- min(1, (1-exp(-d1/N.total))/(1-(exp(-L.0/N.total))))
  R2h.single.max <- data.frame(cbind(R2h.single.max.value, R2h.single.max.lb, R2h.single.max.ub))  # OUT
  colnames(R2h.single.max) <- c("R2h.single.max", "CI lower limit", "CI upper limit")
  rownames(R2h.single.max) <- c(" ")
    
  # corrs both treatment groups
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
    list(Data.Analyze=wide, Obs.Per.Trial=Obs.per.trial, Trial.Spec.Results=Trial.Spec.Results,  
         R2ht=R2ht, R2h.ind=R2h.single, Cor.Endpoints=Cor.Endpoints, Residuals=Residuals, 
         Call=match.call())   
  
  class(fit) <- "MixedContContIT"
  fit  
}
