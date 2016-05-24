FixedBinContIT <- function(Dataset, Surr, True, Treat, Trial.ID, Pat.ID, 
                    Model=c("Full"), Weighted=TRUE, Min.Trial.Size=2, Alpha=.05, Number.Bootstraps=50, 
                    Seed=sample(1:1000, size=1)){
  
  if ((Model==c("Full") | Model==c("Reduced") | Model==c("SemiReduced"))==FALSE) {stop ("The specification of the Model=c(\"...\") argument of the call is incorrect. Use either Model=c(\"Full\"), Model=c(\"Reduced\"), or Model=c(\"SemiReduced\").")}     
  
  distribution.S <- gaussian
  link.choice.S <- 'identity'
  distribution.T <- binomial
  link.choice.T <- 'logit'
    
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
Model.0 <- glm(True ~ -1 + as.factor(Trial.ID) + as.factor(Trial.ID):Treat, family=distribution.T(link=link.choice.T), data=wide)      
Model.1 <- glm(True ~ -1 + as.factor(Trial.ID) + as.factor(Trial.ID):Treat + as.factor(Trial.ID):Surr,  # LS!
               family=distribution.T(link=link.choice.T), data=wide)
L0 <- -2 * logLik(Model.0)[1]    
L1 <- -2 * logLik(Model.1)[1]    

  # R2h.ind (single-cluster based) and CI
g2 <- -(L1-L0)
R2h.ind.value <- 1 - exp(-g2/N.total)
k1 <- qchisq(Alpha, 1, g2)
d1 <- qchisq((1-Alpha), 1, g2)
R2h.ind.lb <- max(0, 1-exp(-k1/N.total)) 
R2h.ind.ub <- min(1, 1-exp(-d1/N.total))
R2h.ind <- data.frame(cbind(R2h.ind.value, R2h.ind.lb, R2h.ind.ub))   # OUT
colnames(R2h.ind) <- c("R2h.ind", "CI lower limit", "CI upper limit")
rownames(R2h.ind) <- c(" ") 

# R2h.ind.max (single-cluster based) and CI
#Model.0 <- glm(True ~ 1, family=distribution.T(link=link.choice.T), data=wide)      
#L.0 <- -2 * logLik(Model.0)[1]    
#R2h.ind.max.value <- R2h.ind.value/(1-(exp(-L.0/N.total)))   
#R2h.ind.max.lb <- max(0, (1-exp(-k1/N.total))/(1-(exp(-L.0/N.total))))
#R2h.ind.max.ub <- min(1, (1-exp(-d1/N.total))/(1-(exp(-L.0/N.total))))
#R2h.ind.max <- data.frame(cbind(R2h.ind.max.value, R2h.ind.max.lb, R2h.ind.max.ub))  # OUT
#colnames(R2h.ind.max) <- c("R2h.ind.max", "CI lower limit", "CI upper limit")
#rownames(R2h.ind.max) <- c(" ")
  
# Alt.
Model.0.a <- glm(True ~ -1 + as.factor(Trial.ID) + as.factor(Trial.ID):Treat, family=distribution.T(link=link.choice.T), data=wide)      
Model.0.b <- glm(Surr ~ -1 + as.factor(Trial.ID) + as.factor(Trial.ID):Treat, family=distribution.S(link=link.choice.S), data=wide)    
L.0.a <- -2 * logLik(Model.0.a)[1] 
L.0.b <- -2 * logLik(Model.0.b)[1] 
g2 <- -(L1-L0)
L.0 <- min(L.0.a, L.0.b)
R2h.ind.max.value <- g2/L.0

options(warn=-1)
R2h.ind.max.value.boot.all <- NULL
for (z in 1: Number.Bootstraps){  
  aantal <- nrow(wide)
  obs <- c(1:aantal)
  set.seed(Seed)
  Indicator <- sample(obs, aantal, replace=TRUE)  
  Seed <- Seed + 1
  Sample.boot.wide <- data.frame(wide[Indicator,])
  
  Model.0.boot <- glm(True ~ -1 + as.factor(Trial.ID) + as.factor(Trial.ID):Treat, 
                      family=distribution.T(link=link.choice.T), data=Sample.boot.wide)      
  Model.1.boot <- glm(True ~ -1 + as.factor(Trial.ID) + as.factor(Trial.ID):Treat + as.factor(Trial.ID):Surr,  # LS!
                      family=distribution.T(link=link.choice.T), data=Sample.boot.wide)
  L0.boot <- -2 * logLik(Model.0.boot)[1]    
  L1.boot <- -2 * logLik(Model.1.boot)[1]    
  
  g2.boot <- -(L1.boot-L0.boot)
  Model.0.a.boot <- glm(True ~ -1 + as.factor(Trial.ID) + as.factor(Trial.ID):Treat, family=distribution.T(link=link.choice.T), data=Sample.boot.wide)      
  Model.0.b.boot <- glm(Surr ~ -1 + as.factor(Trial.ID) + as.factor(Trial.ID):Treat, family=distribution.S(link=link.choice.S), data=Sample.boot.wide)    
  L.0.a.boot <- -2 * logLik(Model.0.a.boot)[1] 
  L.0.b.boot <- -2 * logLik(Model.0.b.boot)[1] 
  g2.boot <- -(L1.boot-L0.boot)
  L.0.boot <- min(L.0.a.boot, L.0.b.boot)
  R2h.ind.max.value.boot <- g2.boot/L.0.boot
  R2h.ind.max.value.boot.all <- 
    cbind(R2h.ind.max.value.boot.all, R2h.ind.max.value.boot)
}
options(warn=0)

Mean.Boot.CI.R2h.ind <- mean(R2h.ind.max.value.boot.all)
Var.Boot.CI.R2h.ind <- var(as.numeric(R2h.ind.max.value.boot.all))
Sort.CI <- sort(R2h.ind.max.value.boot.all)
Boot.CI.R2h.ind.lb <- max(0, Mean.Boot.CI.R2h.ind + qnorm(Alpha/2)*sqrt(Var.Boot.CI.R2h.ind))    
Boot.CI.R2h.ind.ub <- min(1, Mean.Boot.CI.R2h.ind - qnorm(Alpha/2)*sqrt(Var.Boot.CI.R2h.ind))
R2h.ind.max <- data.frame(cbind(R2h.ind.max.value, sqrt(Var.Boot.CI.R2h.ind),
                                Boot.CI.R2h.ind.lb, Boot.CI.R2h.ind.ub))
colnames(R2h.ind.max) <- c("R2b.ind", "Standard Error", "CI lower limit", "CI upper limit")
rownames(R2h.ind.max) <- c(" ")


  # R2h (multiple-cluster based) and CI                   
trialnames <- unique(wide$Trial.ID)
count <- c(0)
sum.tot <- c(0)

options(warn=-1)
for (i in 1:N.trial) {
   {
    Data.temp <- subset(wide, subset=(Trial.ID==trialnames[i]))  # inf likelihoods als <4 observaties voor Model.1
    if (nrow(Data.temp) >= 4) {    
     if (length(unique(Data.temp$True))!=1){          # do not consider clusters with all same T (inf L)
     if (length(unique(Data.temp$Treat))==2){        # not consider clusters with all same Treat
  
     Model.1.temp <- glm(True ~ Treat, family=distribution.T(link=link.choice.T), data=Data.temp)         
     Model.2.temp <- glm(True ~ Treat + Surr, family=distribution.T(link=link.choice.T), data=Data.temp) 
     L1.temp <- -2 * logLik(Model.1.temp)[1]
     L2.temp <- -2 * logLik(Model.2.temp)[1]
     n.i <- nrow(Data.temp)
     if (L1.temp == Inf | L2.temp == Inf) {cat("Problem: infinite likelihood for one of the models of cluster ", i, "\n")}
    
     G.i <- -(L2.temp-L1.temp)
     
     if (G.i != 0){
       count <- count + 1    
       sum.tot.part <- exp(-G.i/n.i)    # G_i niet ^2
       sum.tot <- sum.tot + sum.tot.part}
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
    
     Model.1.temp <- glm(eval(parse(text=command.Model.1.temp)), family=distribution.T(link=link.choice.T), data=Data.temp)         
     Model.2.temp <- glm(eval(parse(text=command.Model.2.temp)), family=distribution.T(link=link.choice.T), data=Data.temp) 
     L1.temp <- -2 * logLik(Model.1.temp)[1]
     L2.temp <- -2 * logLik(Model.2.temp)[1]
     n.i <- nrow(Data.temp)

     if (abs(L1.temp) != Inf | abs(L2.temp) != Inf){   
     G.i <- -(L2.temp-L1.temp)
     if (G.i != 0){
       count <- count + 1    
       sum.tot.part <- exp(-((G.i) / n.i))    
       sum.tot <- sum.tot + sum.tot.part}
     }
     }} }
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

# Per trial R2_h  

R2h.per.trial.lb.vec <- R2h.per.trial.Name.vec <- R2h.per.trial.ub.vec <- R2h.per.trial.hier.vec <- NULL
trialnames <- unique(wide$Trial.ID)

options(warn=-1)
for (i in 1:N.trial) {
  Data.temp <- subset(wide, subset=(Trial.ID==trialnames[i]))  # inf likelihoods als <4 observaties voor Model.1
  if (nrow(Data.temp) >= 4) {    
    if (length(unique(Data.temp$True))!=1){          # do not consider clusters with all same T (inf L)
      if (length(unique(Data.temp$Treat))==2){        # not consider clusters with all same Treat
        
        R2h.per.trial.Name <- trialnames[i]
        # R2h.per.trial 
        Model.0.hier <- glm(True ~ 1 + Treat, 
                            family=distribution.T(link=link.choice.T), data=Data.temp)      
        Model.1.hier <- glm(True ~ 1 + Treat + Surr,  # LS!
                            family=distribution.T(link=link.choice.T), data=Data.temp)
        L0.hier <- -2 * logLik(Model.0.hier)[1]    
        L1.hier <- -2 * logLik(Model.1.hier)[1]    
        N.total.hier <- dim(Data.temp)[1]
        
        g2.hier <- -(L1.hier - L0.hier)
        
        if (g2.hier != 0){
        R2h.per.trial.hier <- 1 - exp(-g2.hier/N.total.hier)
        k1 <- qchisq(Alpha, 1, g2.hier)
        d1 <- qchisq((1-Alpha), 1, g2.hier)
        R2h.per.trial.lb <- max(0, 1-exp(-k1/N.total.hier)) 
        R2h.per.trial.ub <- min(1, 1-exp(-d1/N.total.hier))
        
        R2h.per.trial.hier.vec <- cbind(R2h.per.trial.hier.vec, R2h.per.trial.hier)
        R2h.per.trial.lb.vec <- cbind(R2h.per.trial.lb.vec, R2h.per.trial.lb)
        R2h.per.trial.ub.vec <- cbind(R2h.per.trial.ub.vec, R2h.per.trial.ub)
        R2h.per.trial.Name.vec <- cbind(R2h.per.trial.Name.vec, R2h.per.trial.Name)
      }}}}}
options(warn=0)

R2h.By.Trial <- cbind(t(R2h.per.trial.Name.vec), t(R2h.per.trial.hier.vec), t(R2h.per.trial.lb.vec), 
        t(R2h.per.trial.ub.vec))
colnames(R2h.By.Trial) <- c("TrialID", "R2h", "R2h_low", "R2h_up")
rownames(R2h.By.Trial) <- NULL 

fit <- 
  list(Data.Analyze=wide, Obs.Per.Trial=Obs.per.trial, Trial.Spec.Results=Trial.Spec.Results,  
       R2ht=R2ht, R2h.ind=R2h.ind, R2h=R2h.cluster.based, R2b.ind=R2h.ind.max, 
       R2h.Ind.By.Trial=R2h.By.Trial, Call=match.call()) 

 
  class(fit) <- "FixedBinContIT"
  fit  
  
}





summary.FixedBinBinIT <- summary.FixedBinContIT <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Data summary and descriptives")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\nTotal number of trials: ", nrow(Object$Obs.Per.Trial))
  cat("\nTotal number of patients: ", dim(Object$Data.Analyze)[1])
  cat("\nM(SD) patients per trial: ", format(round(mean((Object$Obs.Per.Trial$Obs.per.trial)), 4), nsmall = 4), " (", format(round(sd((Object$Obs.Per.Trial$Obs.per.trial)), 4), nsmall = 4), ")", 
      "  [min: ", min((Object$Obs.Per.Trial$Obs.per.trial)), "; max: ",  max((Object$Obs.Per.Trial$Obs.per.trial)), "]", sep="")
  cat("\nTotal number of patients in experimental treatment group: ", length(Object$Data.Analyze$Treat[Object$Data.Analyze$Treat==1]), 
      "\nTotal number of patients in control treatment group: ", length(Object$Data.Analyze$Treat[Object$Data.Analyze$Treat!=1])) 
  
  cat("\n\n\n# Information-theoretic surrogacy estimates summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\n")
  cat("Trial-level surrogacy (R2_ht): \n")
  print(format(round(Object$R2ht, 4), nsmall = 4))
#  cat("\nIndividual-level surrogacy (R2_h): \n")
#  print(format(round(Object$R2h, 4), nsmall = 4))
  cat("\nIndividual-level surrogacy (R2_h.ind): \n")
  print(format(round(Object$R2h.ind, 4), nsmall = 4))
  cat("\nIndividual-level surrogacy taking max. bound into account (R2_b.ind): \n")
  print(format(round(Object$R2b.ind, 4), nsmall = 4))
  
}





plot.FixedBinBinIT <- plot.FixedBinContIT <- function(x, Trial.Level=TRUE, Weighted=TRUE, Indiv.Level.By.Trial=TRUE, 
                               Xlab.Indiv, Ylab.Indiv, Xlab.Trial, Ylab.Trial, Main.Trial, Main.Indiv, 
                               Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)),...){
  
  Object <- x 
  if (Trial.Level==TRUE){ 
    if (missing(Xlab.Trial)) {Xlab.Trial <- expression(paste("Treatment effect on the Surrogate endpoint ", (alpha[i])))}
    if (missing(Ylab.Trial)) {Ylab.Trial <- expression(paste("Treatment effect on the True endpoint  ",(beta[i])))}
    if (missing(Main.Trial)) {Main.Trial <- c("Trial-level surrogacy")}
    dev.new()
    par=Par
    
    if (Weighted==TRUE){
      plot(Object$Trial.Spec.Results$Treatment.S, Object$Trial.Spec.Results$Treatment.T, 
           cex=(Object$Trial.Spec.Results$Obs.per.trial/max(Object$Trial.Spec.Results$Obs.per.trial))*8, 
           xlab=Xlab.Trial, ylab=Ylab.Trial, main=Main.Trial,...)
      abline(lm(Object$Trial.Spec.Results$Treatment.T ~ Object$Trial.Spec.Results$Treatment.S))}
    
    if (Weighted==FALSE){
      plot(Object$Trial.Spec.Results$Treatment.S, Object$Trial.Spec.Results$Treatment.T, 
           xlab=Xlab.Trial, ylab=Ylab.Trial, main=Main.Trial, ...)
      abline(lm(Object$Trial.Spec.Results$Treatment.T ~ Object$Trial.Spec.Results$Treatment.S))}
  }
  
  if (Indiv.Level.By.Trial==TRUE){ 
    if (missing(Xlab.Indiv)) {Xlab.Indiv <- expression(paste("R"[h.ind]^{2}))}
    if (missing(Ylab.Indiv)) {Ylab.Indiv <- "Trial"}
    if (missing(Main.Indiv)) {Main.Indiv <- c("Individual-level surrogacy")}
    dev.new()
    par=Par
    plot(y=(1:dim(Object$R2h.Ind.By.Trial)[1]), x=Object$R2h.Ind.By.Trial[,2], yaxt = "n",
         xlab=Xlab.Indiv, ylab=Ylab.Indiv, main=Main.Indiv, xlim=c(0, 1), ...)
    axis(2, at=(1:dim(Object$R2h.Ind.By.Trial)[1]), labels=Object$R2h.Ind.By.Trial[,1])
    for (i in 1: dim(Object$R2h.Ind.By.Trial)[1]){
      segments(y0 = i, x0 = Object$R2h.Ind.By.Trial[i,3], y1 = i, x1 = Object$R2h.Ind.By.Trial[i,4], lty=2)
    }
    
  }    
} 
