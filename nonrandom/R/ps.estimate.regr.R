ps.estimate.regr <- function(data,
                             name.data.regr,
                             regr,
                             name.regr,
                             treat,
                             name.treat,
                             resp,
                             name.resp,                         
                             lr.form,
                             family,
                             ...)
{

  ## ####################
  ## fit regression model

  ## conditional odds ratio
  y.model <- glm(lr.form,
                 family = family,
                 data = data,
                 ...)
  
  effect.lr <-
    summary(y.model)$coeff[rownames(summary(y.model)$coeff) == name.treat,1]
  
  se.effect.lr  <-
    summary(y.model)$coeff[rownames(summary(y.model)$coeff) == name.treat,2]
  
  
  ## marginal odds ratio      
  if (family == "binomial"){
    
    effect.lr <-
      exp(summary(y.model)$coeff[rownames(summary(y.model)$coeff) == name.treat,1])
    
    
    p1.lr <- 1/(1+exp(-(as.vector(summary(y.model)$coeff[,1])%*%
                        t(as.data.frame(cbind(rep(1, dim(data)[1]),
                                              rep(1, dim(data)[1]), regr))))))
    p0.lr <- 1/(1+exp(-(as.vector(summary(y.model)$coeff[,1])%*%
                        t(as.data.frame(cbind(rep(1, dim(data)[1]),
                                              rep(0, dim(data)[1]), regr))))))
    
    effect.lr.marg <-  (mean(p1.lr, na.rm=TRUE)/(1-mean(p1.lr, na.rm=TRUE))) /
      (mean(p0.lr, na.rm=TRUE)/(1-mean(p0.lr, na.rm=TRUE)))
    
    
    ## variance of marginal odds ratio
    logit.p1.lr <- mean(p1.lr, na.rm=TRUE)*(1-mean(p1.lr, na.rm=TRUE))
    logit.p0.lr <- mean(p0.lr, na.rm=TRUE)*(1-mean(p0.lr, na.rm=TRUE))
    
    weight1.lr <- (p1.lr*(1-p1.lr))[!is.na(p1.lr)]
    weight0.lr <- (p0.lr*(1-p0.lr))[!is.na(p0.lr)]
    
    x1 <- as.matrix(cbind(rep(1, dim(data)[1]),
                          rep(1, dim(data)[1]),
                          regr))[!is.na(p1.lr),]
    x0 <- as.matrix(cbind(rep(1, dim(data)[1]),
                          rep(0, dim(data)[1]),
                          regr))[!is.na(p0.lr),]
    
    der1 <- apply(x1, 2, function(x) x*weight1.lr)
    der0 <- apply(x0, 2, function(x) x*weight0.lr)
    
    cov.pi1 <- der1 %*% vcov(y.model) %*% t(der1)
    cov.pi0 <- der0 %*% vcov(y.model) %*% t(der0) 
    
    se.lr.marg <-
      sqrt((1/logit.p1.lr^2)*(sum(cov.pi1)/dim(x1)[1]^2) + (1/logit.p0.lr^2)*(sum(cov.pi0)/dim(x0)[1]^2) -
           2*(1/logit.p1.lr)*(1/logit.p0.lr)*(cov(p0.lr[1,!is.na(p0.lr)],p1.lr[1,!is.na(p1.lr)])/
                                              dim(data)[1])) ## geaendert am 22.10. : dim(data)[1]^2
    
  }else{
    
    effect.lr.marg <- effect.lr
    se.lr.marg <- se.effect.lr
    
  }

  lr.estimation <- list(effect       = effect.lr,
                        effect.marg  = effect.lr.marg,
                        se           = se.effect.lr,
                        se.marg      = se.lr.marg,
                        regr.formula = lr.form,
                        regr.model   = y.model)

  return(lr.estimation)
  
}

