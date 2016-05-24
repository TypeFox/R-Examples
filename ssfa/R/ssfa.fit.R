ssfa.fit <- function(y, x, w, intercept = TRUE, pars = NULL, sigmau2_sar =NULL, 
                     par_rho=TRUE, form = "cost", fun = "hnormal"){
  options(warn=-1)
  p <- ncol(x)
  X <- as.matrix(x)
    
  #neighbours matrix
  w <- as.matrix(w)
  list_w <- mat2listw(w)

  
  # type of form
  if (form == "cost") sc = -1
  if (form == "production") sc = 1

  # intercept
  if (intercept) {
    X <- cbind(1, X)
  }

  data <- data.frame(y, X)
  
  #ols
  ols <- lm(y ~ ., data = data)
  if (intercept) {
  ols <- lm(y ~ .- 1, data = data)
  }
  
  #########starting values############
  
    m2 <- mean(ols$residuals^2)
  
  #3rd moment test
  
    m3 <- mean(ols$residuals^3)
    m3t<- m3/sqrt(6*m2^3/length(y))

  #correct 3rd moment to be negative
  
  if(sc*sum(ols$residuals^3)<0) { 
    m3<-m3 
  } else
  {
    m3<- -0.0001*sc
  }
    
    #su <- var(ols$residuals)*(1-2/pi)
    su2 <- (sc*m3/(sqrt(2/pi)*(1-4/pi)))^(2/3)
    su <- sqrt(su2)
    #sv <- var(ols$residuals)
    sv2 <- m2 - (1-2/pi)*su2
    if(sv2>0) {
      sv2<-sv2
  }  else
  {
    sv2 <- 0.0001 
  } 
   sv <- sqrt(sv2) 
 
  ### b0
  if (intercept) {
  b0 <- as.matrix(ols$coefficients[1]+sc*(sqrt(2/pi))*su)
  } 
  
    
  # spatial dependence
  
  if(par_rho == TRUE)
  {
  sarlm <- errorsarlm(y ~ . , data=data, listw=list_w, method="eigen", quiet=NULL)
  sigmau2_sar <- var(ols$residuals-sarlm$residuals)	
  if (is.null(pars)) {
    pars <- c(sarlm$lambda, su2, sv2, ols$coefficients)
    if (intercept) {
    pars <- c(sarlm$lambda, su2, sv2, b0, ols$coefficients[2:length(ols$coefficients)])
      names(pars) <- c("rho", "sigmau2_dmu", "sigmav2", "Intercept", names(pars)[5:length(pars)])
      
      #type of inefficiency distribution
      maxlik = L_hNV_rho
      
      mod <- maxNR(L_hNV_rho, start=pars, X=X, y=y, w=w, sigmau2_sar=sigmau2_sar, sc=sc, 
                   finalHessian=T)
    
      coef <- mod$estimate
      
      sigma2 <- coef["sigmau2_dmu"] + sigmau2_sar +  coef["sigmav2"]
      names(sigma2) <- "sigma2"
       
      if(coef["sigmau2_dmu"]<=0)
      {
        LogLik <- logLik(ols)
      }
      else
      {
        LogLik <-  LogLik <- mod$maximum
      }
        
      ret <- list(y = y, x = x, X=X, w=w, coef = coef, sc = sc, hess = mod$hessian, logLik = LogLik,
                  ols = ols, sigmau2_dmu=coef["sigmau2_dmu"], sigmav2=coef["sigmav2"], 
                  sigmau2_sar=sigmau2_sar, sigma2=sigma2, fun=fun, rho=coef["rho"], list_w=list_w)
    }
   }
  }
  if(par_rho == FALSE) ### return SFA
  {
    rho <- 0
    sigmau2_sar<-0
    if (is.null(pars)) {
    pars <- c(su2, sv2, ols$coefficients)
    if (intercept) {
      pars <- c(su2, sv2, b0, ols$coefficients[2:length(ols$coefficients)])
      names(pars) <- c("sigmau2", "sigmav2", "Intercept", names(pars)[4:length(pars)])
      
      #type of inefficiency distribution
      maxlik = L_hNV
      
      mod <- maxNR(L_hNV, start=pars, X=X, y=y, sc=sc, finalHessian=T)
      
        
      coef <- mod$estimate

      sigma2 <- coef["sigmau2"] + coef["sigmav2"]
      names(sigma2) <- "sigma2"
          
      LogLik <- mod$maximum
      
       ret <- list(y = y, x = x, X=X, w=w, coef = coef, sc = sc, hess = mod$hessian, logLik = LogLik,
                   ols = ols, sigmau2=coef["sigmau2"], sigmav2=coef["sigmav2"], 
                   sigmau2_sar=sigmau2_sar, sigma2=sigma2, fun=fun, rho=rho, list_w=list_w)

    }
   }
  }
    
  class(ret) <- "ssfa"
  return(ret)
}


