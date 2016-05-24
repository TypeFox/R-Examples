summary.ssfa <- function(object,...){
  coef <- object$coef
  var_beta <- diag(abs(solve(object$hess)))  
  zvalue <- round(coef/sqrt(var_beta), digits=6)
  pvalue_coeff <- round(2*pnorm(abs(zvalue), lower.tail = FALSE), digits=6) 
  
  if(object$rho==0)
  {
    coef <-c(coef["Intercept"], coef[names(coef)[4:length(coef)]], coef["sigmau2"], coef["sigmav2"])
    sqrt_var_beta <-sqrt(c(var_beta["Intercept"], var_beta[names(object$coef)[4:length(object$coef)]], var_beta["sigmau2"], var_beta["sigmav2"]))
    zvalue <-c(zvalue["Intercept"], zvalue[names(object$coef)[4:length(object$coef)]], zvalue["sigmau2"], zvalue["sigmav2"])
    pvalue_coeff <-c(pvalue_coeff["Intercept"], pvalue_coeff[names(object$coef)[4:length(object$coef)]], pvalue_coeff["sigmau2"], pvalue_coeff["sigmav2"])
    
    coef.table <- cbind(coef, sqrt_var_beta, zvalue, pvalue_coeff)
    colnames(coef.table) <- c("Estimate", "Std. Error","z value","Pr(>|z|)")
    lambda <- sqrt(coef["sigmau2"]/coef["sigmav2"])
    
    l_ssfa <- object$logLik
    l_ols <- logLik(object$ols)
    
    L <- round(2*(l_ssfa-l_ols)[1], digits = 3)
    names(L) <- "LR-Test"
    pvalue_eff <-  pchisq(abs(L), 1,lower.tail = FALSE) / 2  
    
    names(pvalue_eff) <- "p-value"
    test.table <- rbind(l_ssfa, l_ols)
    dimnames(test.table) <- list(NULL, "Value Log-Lik")
    row.names(test.table) <- c("ssfa", "ols")
    wert.table <- list(L=L, pvalue_eff=pvalue_eff)
    rho <- round(object$rho, digits=6)
    sigmau2_sar <- object$sigmau2_sar
    mean.eff <- mean(eff.ssfa(object))
    names(mean.eff) <- "mean efficiency"
    moran_test <- moran.test(residuals.ssfa(object), object$list_w)
    
    AIC_ssfa <- -2*l_ssfa +  2*length(object$coef)
    AIC_ols  <- -2*l_ols  +  2*length(object$ols$coef)
    
    sigma2t <- sum(residuals.ssfa(object)^2)/length(residuals.ssfa(object))
    
    ret <- list(coef.table=coef.table, test.table = test.table, pvalue_coeff=pvalue_coeff, pvalue_eff = 
                  pvalue_eff, df = df, LR = L, mean.eff = mean.eff, rho=rho, moran_test=moran_test, 
                AIC_ssfa=AIC_ssfa, AIC_ols=AIC_ols, lambda=lambda, sigma2t=sigma2t, sigmau2_sar=sigmau2_sar)
    
  }
  
  if(object$rho!=0)
  {
    coef <-c(coef["Intercept"], coef[names(coef)[5:length(coef)]], coef["sigmau2_dmu"], coef["sigmav2"])
    sqrt_var_beta <-sqrt(c(var_beta["Intercept"], var_beta[names(object$coef)[5:length(object$coef)]], var_beta["sigmau2_dmu"], var_beta["sigmav2"]))
    zvalue <-c(zvalue["Intercept"], zvalue[names(object$coef)[5:length(object$coef)]], zvalue["sigmau2_dmu"], zvalue["sigmav2"])
    pvalue_coeff <-c(pvalue_coeff["Intercept"], pvalue_coeff[names(object$coef)[5:length(object$coef)]], pvalue_coeff["sigmau2_dmu"], pvalue_coeff["sigmav2"])
    
    coef.table <- cbind(coef, sqrt_var_beta, zvalue, pvalue_coeff)
    colnames(coef.table) <- c("Estimate", "Std. Error","z value","Pr(>|z|)")
    if(coef["sigmau2_dmu"]<=0)
    {
      lambda <- 0
    }
    else
    {
      lambda <- coef["sigmau2_dmu"]/coef["sigmav2"]
    }
  
    l_ssfa <- object$logLik
    l_ols <- logLik(object$ols)
    L <- round(2*(l_ssfa-l_ols)[1], digits = 3)
    names(L) <- "LR-Test"
    pvalue_eff <-  pchisq(abs(L), 1,lower.tail = FALSE) / 2  
    names(pvalue_eff) <- "p-value"
    test.table <- rbind(l_ssfa, l_ols)
    dimnames(test.table) <- list(NULL, "Value Log-Lik")
    row.names(test.table) <- c("ssfa", "ols")
    wert.table <- list(L=L, pvalue_eff=pvalue_eff)
    rho <- round(object$rho, digits=6)
    sigmau2_sar <- object$sigmau2_sar
    mean.eff <- mean(eff.ssfa(object))
    names(mean.eff) <- "mean efficiency"
    moran_test <- moran.test(residuals.ssfa(object), object$list_w)
    
    AIC_ssfa <- -2*l_ssfa + 2*2*length(object$coef)
    AIC_ols  <- -2*l_ols + 2*length(object$ols$coef)
    
    sigma2t <- sum(residuals.ssfa(object)^2)/length(residuals.ssfa(object))
    
    sigmau2 <- object$sigmau2_sar + object$sigmau2_dmu
    
    ret <- list(coef.table=coef.table, test.table = test.table, pvalue_coeff=pvalue_coeff, pvalue_eff = 
                  pvalue_eff, df = df, LR = L, mean.eff = mean.eff, rho=rho, moran_test=moran_test, 
                AIC_ssfa=AIC_ssfa, AIC_ols=AIC_ols, lambda=lambda, sigma2t=sigma2t, sigmau2_sar=sigmau2_sar,
                sigmau2=sigmau2)
    
  }
  
  ### class
  class(ret) <- "summary.ssfa"
  return(ret)
}
