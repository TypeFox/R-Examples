Single.Trial.RE.AA <- function(Dataset, Surr, True, Treat, Pat.ID, Alpha=.05, Number.Bootstraps=500, Seed=sample(1:1000, size=1)){
  
  Surr <- Dataset[,paste(substitute(Surr))]
  True <- Dataset[,paste(substitute(True))]
  Treat <- Dataset[,paste(substitute(Treat))]
  Trial.ID <- rep(1, nrow(Dataset))
  Pat.ID <- Dataset[,paste(substitute(Pat.ID))]
  
  Data.Proc <- .Data.Processing(Dataset=Dataset, Surr=Surr, True=True, Treat=Treat, Trial.ID=Trial.ID, Pat.ID=Pat.ID, Min.Trial.Size=0)
  wide <- Data.Proc$wide
  dataS <- Data.Proc$dataS
  dataT <- Data.Proc$dataT
  Data.analyze <- Data.Proc$Data.analyze
  N.total <- Data.Proc$N.total

  model12 <- lm(cbind(wide$Surr, wide$True)~ wide$Treat, data=wide) 
  Residuals <- data.frame(model12$residuals)
  colnames(Residuals) <- c("Surr", "True")
  alpha <- model12$coefficients[2,1]
  alpha_se <- summary(model12)$"Response wide$Surr"$coefficients[2,2]
  alpha_lb <- alpha - (qt(c(1-Alpha/2), df=N.total, lower.tail=TRUE)*alpha_se)
  alpha_ub <- alpha + (qt(c(1-Alpha/2), df=N.total, lower.tail=TRUE)*alpha_se)
  alpha_results <- data.frame(cbind(alpha, alpha_se , alpha_lb, alpha_ub))
  colnames(alpha_results) <- c("Alpha", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(alpha_results) <- c(" ")
  
  beta <- model12$coefficients[2,2]
  beta_se <- summary(model12)$"Response wide$True"$coefficients[2,2]
  beta_lb <- beta - (qt(c(1-Alpha/2), df=N.total, lower.tail=TRUE)*beta_se)
  beta_ub <- beta + (qt(c(1-Alpha/2), df=N.total, lower.tail=TRUE)*beta_se)
  beta_results <- data.frame(cbind(beta, beta_se , beta_lb, beta_ub))
  colnames(beta_results) <- c("Beta", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(beta_results) <- c(" ")
  
  # CI AA
  rho_z <- var(model12$residuals[,2], model12$residuals[,1])/ sqrt(var(model12$residuals[,2])*var(model12$residuals[,1]))
  Z <- .5*log((1+rho_z)/(1-rho_z)) 
  rho_lb <- max(0, (exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))-1)/(exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))+1))
  rho_ub <- min(1, (exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))-1)/(exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))+1))
  rho_sd <- sqrt((1-rho_z**2)/(N.total-2))
  rho_results_FishZ <- data.frame(cbind(rho_z, rho_sd , rho_lb, rho_ub))
  colnames(rho_results_FishZ) <- c("AA (gamma)", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(rho_results_FishZ) <- c(" ")
  
  # RE delta
  RE <- beta/alpha
  cov_a_b_s <- vcov(model12)[4,2]  
  var_alpha <- vcov(model12)[2,2]  
  var_beta <- vcov(model12)[4,4]   
  estm <- as.numeric(cbind(beta, alpha))
  estv <- matrix(c(var_beta, cov_a_b_s, cov_a_b_s, var_alpha), nrow=2)
  std <- deltamethod (~ (x1/x2), estm, estv)   
  RE_low <- as.numeric(RE+(qnorm(Alpha/2)*std))
  RE_high <- as.numeric(RE+(qnorm(1-Alpha/2)*std))
  RE_results_Delta <- data.frame(cbind(RE, std, RE_low, RE_high))   
  colnames(RE_results_Delta) <- c("RE", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(RE_results_Delta) <- c(" ")  
  
  # RE Fieler 
  cov_a_b_s <- vcov(model12)[4,2]  
  var_alpha <- vcov(model12)[2,2]  
  var_beta <- vcov(model12)[4,4]   
  A <- (alpha * beta) - (((qnorm(1-Alpha/2))^2)*cov_a_b_s)
  B <- (alpha^2) - (((qnorm(1-Alpha/2))^2)*var_alpha)
  C <- (beta^2) - (((qnorm(1-Alpha/2))^2)*var_beta)
  
  if ((A^2)-(B*C) < 0) {
  RE_results_Fiel <- data.frame(cbind(RE, std, NaN, NaN))   
  colnames(RE_results_Fiel) <- c("RE", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(RE_results_Fiel) <- c(" ")}  
  
  if ((A^2)-(B*C) > 0) {
  try(RE_low <- ((A - sqrt((A^2)-(B*C)))/B), silent=TRUE)       
  try(RE_high <- ((A + sqrt((A^2)-(B*C)))/B), silent=TRUE)
  var <- ((beta/alpha)^2) * ((var_beta / (beta^2)) + (var_alpha / (alpha^2)) - (2 * (cov_a_b_s/(alpha*beta)))) 
  std <- sqrt(var)
  RE_results_Fiel <- data.frame(cbind(RE, std, RE_low, RE_high))   
  colnames(RE_results_Fiel) <- c("RE", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(RE_results_Fiel) <- c(" ") 
  }
  
  d.size <- dim(wide)[1] 
  obs <- c(1:d.size)   
  k <- Number.Bootstraps  
  RE_boot <- rho_z_boot <- alpha_boot <- beta_boot <- as.vector(NULL)
    for (i in 1:k){
    set.seed(Seed+i)
    index <- sample(obs, d.size, replace=TRUE)
    sample <- data.frame(wide[index,])                             
    sample <- na.exclude(sample[order(sample$Pat.ID),])            
    model1 <- lm(cbind(sample$Surr, sample$True)~ sample$Treat, data=sample)
    alpha_boot[i] <- model1$coefficients[2,1]
    beta_boot[i] <- model1$coefficients[2,2]
    rho_z_boot[i] <- var(model1$residuals[,2], model1$residuals[,1])/ sqrt(var(model1$residuals[,2])*var(model1$residuals[,1]))
    RE_boot[i] <- beta_boot[i] / alpha_boot[i]
  }
 
  res_RE <- studres(lm(RE_boot ~ 1))
  if ((max(abs(res_RE)) > qt(c(1-(Alpha/(2*Number.Bootstraps))), df=(Number.Bootstraps-1-1), lower.tail=TRUE))==TRUE){
  cat("\nWarning: There were outliers in the bootstrapped RE sample.")
  cat("\nThe bootstrap-based standard error and/or confidence interval of RE may not be")
  cat("\ntrustworthy. The following observations (in the sample of bootstrapped RE values)")
  cat("\nare outliers (using abs(", qt(c(1-(Alpha/(2*Number.Bootstraps))), df=(Number.Bootstraps-1-1), lower.tail=TRUE), ") as the critical value):\n", sep="")
  print(res_RE[abs(res_RE) > qt(c(1-(Alpha/(2*Number.Bootstraps))), df=(Number.Bootstraps-1-1), lower.tail=TRUE)])
  }
  
  RE_CIs <- quantile(RE_boot, probs=c(Alpha/2, 1-Alpha/2))                                                                     
  RE_results_Boot <- data.frame(cbind(RE, sd(RE_boot), RE_CIs[1], RE_CIs[2]))   
  colnames(RE_results_Boot) <- c("RE", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(RE_results_Boot) <- c(" ")
  
  rho_CIs <- quantile(rho_z_boot, probs=c(Alpha/2, 1-Alpha/2)) 
  rho_results_Boot <- data.frame(cbind(rho_z, sd(rho_z_boot), rho_CIs[1], rho_CIs[2]))   
  colnames(rho_results_Boot) <- c("AA (gamma)", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(rho_results_Boot) <- c(" ")

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
  list(Data.Analyze=wide, Alpha=alpha_results, Beta=beta_results, RE.Delta=RE_results_Delta, RE.Fieller=RE_results_Fiel, RE.Boot=RE_results_Boot, RE.Boot.Samples=RE_boot, AA=rho_results_FishZ, AA.Boot=rho_results_Boot, AA.Boot.Samples=rho_z_boot,
       Cor.Endpoints=Cor.Endpoints, Residuals=Residuals, Call=match.call())  

class(fit) <- "Single.Trial.RE.AA"
fit

}