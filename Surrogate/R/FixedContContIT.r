FixedContContIT <- function(Dataset, Surr, True, Treat, Trial.ID, Pat.ID, 
                    Model=c("Full"), Weighted=TRUE, Min.Trial.Size=2, Alpha=.05, Number.Bootstraps=500, 
                    Seed=sample(1:1000, size=1)){
  
  if ((Model==c("Full") | Model==c("Reduced") | Model==c("SemiReduced"))==FALSE) {stop ("The specification of the Model=c(\"...\") argument of the call is incorrect. Use either Model=c(\"Full\"), Model=c(\"Reduced\"), or Model=c(\"SemiReduced\").")}     
  
  distribution.S <- gaussian
  link.choice.S <- 'identity'
  distribution.T <- gaussian
  link.choice.T <- 'identity'
    
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

# R2_ht          

  if (Model==c("Full")|Model==c("SemiReduced")){
    Model.S <- glm(outcome ~ -1 + as.factor(Trial.ID) + as.factor(Trial.ID):Treat, data=dataS, family=distribution.S(link=link.choice.S)) 
    Model.T <- glm(outcome ~ -1 + as.factor(Trial.ID) + as.factor(Trial.ID):Treat, data=dataT, family=distribution.T(link=link.choice.T))
  }
  
  if (Model==c("Reduced")){
    Model.S <- glm(outcome ~ as.factor(Trial.ID):Treat, data=dataS, family=distribution.S(link=link.choice.S)) 
    Model.T <- glm(outcome ~ as.factor(Trial.ID):Treat, data=dataT, family=distribution.T(link=link.choice.T))
  }  
  
 
if (Model==c("Full")| Model==c("SemiReduced")){  
Intercept.S <- coef(Model.S)[1:N.trial]         
Treatment.S <- coef(Model.S)[(1+N.trial):(2*N.trial)] 
Intercept.T <- coef(Model.T)[1:N.trial]
Treatment.T <- coef(Model.T)[(1+N.trial):(2*N.trial)] 
Trial.Spec.Results <- data.frame(Obs.per.trial$Trial, Obs.per.trial$Obs.per.trial, Intercept.S, Intercept.T, Treatment.S, Treatment.T)
colnames(Trial.Spec.Results) <- c(NULL, "Trial", "Obs.per.trial", "Intercept.S", "Intercept.T", "Treatment.S", "Treatment.T")
rownames(Trial.Spec.Results) <- NULL 
}

if (Model==c("Reduced")){  
Treatment.S <- coef(Model.S)[2:(N.trial+1)]
Treatment.T <- coef(Model.T)[2:(N.trial+1)]
Trial.Spec.Results <- data.frame(Obs.per.trial$Trial, Obs.per.trial$Obs.per.trial, Treatment.S, Treatment.T)
colnames(Trial.Spec.Results) <- c(NULL, "Trial", "Obs.per.trial", "Treatment.S", "Treatment.T")
rownames(Trial.Spec.Results) <- NULL
}

# residuals
Residuals.Model.S <- Model.S$residuals
Residuals.Model.T <- Model.T$residuals
Residuals <- data.frame(Residuals.Model.T, Residuals.Model.S) 
rownames(Residuals) <- NULL
Residuals <- cbind(wide$Pat.ID, Residuals)
colnames(Residuals) <- c("Pat.ID", "Residuals.Model.T", "Residuals.Model.S")    
  
  
if (Model==c("Full")){                                       
if (Weighted==FALSE){
Model.1 <- glm(Treatment.T ~ Intercept.S + Treatment.S, family=gaussian)  
L1 <- -2 * logLik(Model.1)[1]
 }
if (Weighted==TRUE){
Model.1 <- glm(Treatment.T ~ Intercept.S + Treatment.S, family=gaussian, weights=Obs.per.trial$Obs.per.trial)  
L1 <- -2 * logLik(Model.1)[1]
 } 
    }
if (Model==c("Reduced")| Model==c("SemiReduced")){
if (Weighted==FALSE){
Model.1 <- glm(Treatment.T ~ Treatment.S, family=gaussian)  
L1 <- -2 * logLik(Model.1)[1]
 }
if (Weighted==TRUE){
Model.1 <- glm(Treatment.T ~ Treatment.S, family=gaussian, weights=Obs.per.trial$Obs.per.trial)  
L1 <- -2 * logLik(Model.1)[1]
 } 
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
  

  # Individual-level surrogacy   
Model.0 <- glm(True ~ -1 + as.factor(Trial.ID) + as.factor(Trial.ID):Treat, 
               family=distribution.T(link=link.choice.T), data=wide)      
Model.1 <- glm(True ~ -1 + as.factor(Trial.ID) + as.factor(Trial.ID):Treat + as.factor(Trial.ID):Surr,  # LS!
               family=distribution.T(link=link.choice.T), data=wide)
L0 <- -2 * logLik(Model.0)[1]    
L1 <- -2 * logLik(Model.1)[1]    

  # R2h.single (single-cluster based) and CI
g2 <- -(L1-L0)
R2h.single.value <- 1 - exp(-g2/N.total)
k1 <- qchisq(Alpha, 1, g2)
d1 <- qchisq((1-Alpha), 1, g2)
R2h.single.lb <- max(0, 1-exp(-k1/N.total)) 
R2h.single.ub <- min(1, 1-exp(-d1/N.total))
R2h.single <- data.frame(cbind(R2h.single.value, R2h.single.lb, R2h.single.ub))   # OUT
colnames(R2h.single) <- c("R2h.ind", "CI lower limit", "CI upper limit")
rownames(R2h.single) <- c(" ") 

# R2h.single.max (single-cluster based) and CI
Model.0 <- glm(True ~ 1, family=distribution.T(link=link.choice.T), data=wide)      
L.0 <- -2 * logLik(Model.0)[1]    
R2h.single.max.value <- R2h.single.value/(1-(exp(-L.0/N.total)))   
R2h.single.max.lb <- max(0, (1-exp(-k1/N.total))/(1-(exp(-L.0/N.total))))
R2h.single.max.ub <- min(1, (1-exp(-d1/N.total))/(1-(exp(-L.0/N.total))))
R2h.single.max <- data.frame(cbind(R2h.single.max.value, R2h.single.max.lb, R2h.single.max.ub))  # OUT
colnames(R2h.single.max) <- c("R2h.single.max", "CI lower limit", "CI upper limit")
rownames(R2h.single.max) <- c(" ")
  
  
  # R2h (multiple-cluster based) and CI                   
trialnames <- unique(wide$Trial.ID)
count <- c(0)
sum.tot <- c(0)
for (i in 1:N.trial) {
   {
    Data.temp <- subset(wide, subset=(Trial.ID==trialnames[i]))  # inf likelihoods als <4 observaties voor Model.1
    if (nrow(Data.temp) >= 4) {    
     if (length(unique(Data.temp$True))!=1){          # do not consider clusters with all same T (inf L)
     if (length(unique(Data.temp$Treat))==2){        # not consider clusters with all same Treat
  
     Model.1.temp <- lm(True ~ Treat, data=Data.temp)         
     Model.2.temp <- lm(True ~ Treat + Surr, data=Data.temp) 
     L1.temp <- -2 * logLik(Model.1.temp)[1]
     L2.temp <- -2 * logLik(Model.2.temp)[1]
     n.i <- nrow(Data.temp)
     if (L1.temp == Inf | L2.temp == Inf) {cat("Problem: infinite likelihood for one of the models of cluster ", i, "\n")}
    
     G.i <- -(L2.temp-L1.temp)
     count <- count + 1    
     sum.tot.part <- exp(-G.i/n.i)    # G_i niet ^2
     sum.tot <- sum.tot + sum.tot.part
       }
     R2h.ind.c.b.val <- 1-(1/count)*(sum.tot) 
            } }}
     }
   
  
    # Bootstrap CI                       
    k <- Number.Bootstraps
    Boot.CI <-rep(0, k)
    for (j in 1:k){
    obs <- c(1:N.total)
    set.seed(Seed)
    Indicator <- sample(obs, N.total, replace=TRUE)  #sample with replacement from 1...Ntotal
    Sample.Boot <- data.frame(wide[Indicator,])  
    Seed <- Seed + 1
    sum.tot <- count <- c(0)    
 
    for (i in 1:N.trial) {
    Data.temp <- subset(Sample.Boot, subset=(Trial.ID==trialnames[i]))  
    if (nrow(Data.temp) >= 4) {                             # inf likelihoods als <4 observaties voor Model.1          
      if (length(unique(Data.temp$True))!=1){          # do not consider clusters with all same T (inf L)
      if (length(unique(Data.temp$Treat))==2){        # not consider clusters with all same Treat      
      command.Model.1.temp <- c("True ~ Treat")
      command.Model.2.temp <- c("True ~ Treat + Surr")
    
     Model.1.temp <- lm(eval(parse(text=command.Model.1.temp)), data=Data.temp)         
     Model.2.temp <- lm(eval(parse(text=command.Model.2.temp)), data=Data.temp) 
     L1.temp <- -2 * logLik(Model.1.temp)[1]
     L2.temp <- -2 * logLik(Model.2.temp)[1]
     n.i <- nrow(Data.temp)

     if (abs(L1.temp) != Inf | abs(L2.temp) != Inf){   
       if (G.i != 0){
         count <- count + 1    
         sum.tot.part <- exp(-((G.i) / n.i))    
         sum.tot <- sum.tot + sum.tot.part}
     }}}}
    R2h.ind.c.b <- 1-(1/count)*(sum.tot)
     Boot.CI[j] <- R2h.ind.c.b
       }
     }
     
     # nonparametric bootstrap normal confidence interval
     Mean.Boot.CI <- mean(Boot.CI)
     Var.Boot.CI <- var(Boot.CI)
     Sort.CI <- sort(Boot.CI)
     Boot.CI.R2.lb <- max(0, R2h.ind.c.b.val + qnorm(Alpha/2)*sqrt(Var.Boot.CI))    
     Boot.CI.R2.ub <- min(1, R2h.ind.c.b.val - qnorm(Alpha/2)*sqrt(Var.Boot.CI))
     R2h.cluster.based <- data.frame(cbind(R2h.ind.c.b.val, sqrt(Var.Boot.CI), Boot.CI.R2.lb, Boot.CI.R2.ub))
     colnames(R2h.cluster.based) <- c("R2h", "Standard Error", "CI lower limit", "CI upper limit")
     rownames(R2h.cluster.based) <- c(" ")
  
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
         R2ht=R2ht, R2h.ind.clust=R2h.cluster.based, R2h.ind=R2h.single, Boot.CI=Boot.CI, Cor.Endpoints=Cor.Endpoints, Residuals=Residuals, Call=match.call())   
  
  class(fit) <- "FixedContContIT"
  fit  
  
}
