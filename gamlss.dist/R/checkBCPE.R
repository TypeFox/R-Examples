checkBCPE <- function(obj=NULL, mu=10, sigma =.1, nu = .5, tau = 2, ...)
 { 
 # ----- function 1
  fun <- function(x, sigma, nu, tau) 
   {
   log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
       c <- exp(log.c) 
     val <- nu*sigma*(x^tau)-x^(tau-1)+(2*sigma*c^2*(1-nu))/tau
     val
   }
  #----- end function 1
  # when the object do not exist
if (is.null(obj))
  { 
     if (tau <= 1) stop(paste("This check does not apply for tau<=1 in BCPE" , "\n"))
     if (nu <= 0 | nu >= 1) stop(paste("There is no minimum turning point for nu<=0 or nu>=1 in BCPE" , "\n"))
     log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
         c <- exp(log.c)
    sigma0 <- ((tau-1)/(c*nu*tau))*((nu*tau)/(2*(1-nu)*(tau-1)))^(1/tau)   
     if (sigma >= sigma0) stop(paste("There is no minimum turning point for sigma equal to ", sigma, " in BCPE" , "\n"))
     lower <- 1/(2*nu*sigma)
     upper <- 1/(nu*sigma) 
     sol<-uniroot(function(x) fun(x, sigma = sigma, nu = nu, tau = tau), low = lower, up = upper, tol = 0.0001)
     z <- -sol$root
     if(length(nu)>1)  y <- ifelse(nu != 0,mu*(nu*sigma*z+1)^(1/nu),mu*exp(sigma*z))
     else   if (nu != 0) y <- mu*(nu*sigma*z+1)^(1/nu) else y <- mu*exp(sigma*z) 
     p <- pBCPE(y, mu=mu, sigma=sigma, nu=nu, tau=tau)
     p
  }
else
  {  # when the object is gamlss
     # chech family == BCPE
      #if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
      if (!("BCPE"%in%obj$family[1])) stop(paste("This check applies only to the BCPE distribution", "\n", ""))
         mu <- fitted(obj)
      sigma <- fitted(obj,"sigma")
         nu <- fitted(obj,"nu")
        tau <- fitted(obj,"tau") 
          N <- length(fitted(obj))
          P <- rep(0,N)
          Y <- rep(0,N)
     Yunder <- rep(TRUE,N)
    for (i in 1:N)
    {
    
      if (tau[i] <= 1 | (nu[i] <= 0 | nu[i] >= 1) )
         { 
          z <- NA
         }
     else
         { 
          log.c <- 0.5*(-(2/tau[i])*log(2)+lgamma(1/tau[i])-lgamma(3/tau[i]))
              c <- exp(log.c)
         sigma0 <- ((tau[i]-1)/(c*nu[i]*tau[i]))*((nu[i]*tau[i])/(2*(1-nu[i])*(tau[i]-1)))^(1/tau[i]) 
          lower <- 1/(2*nu[i]*sigma[i])
          upper <- 1/(nu[i]*sigma[i])
         if ( sigma[i] >= sigma0 ) 
            {
            z <- NA 
            }
         else
            {
             sol<-uniroot(function(x) fun(x, sigma = sigma[i], nu = nu[i], tau = tau[i]), 
             low = lower, up = upper, tol = 0.0001)
             z <- -sol$root
            }
         }
       if (nu[i] != 0) y <- mu[i]*(nu[i]*sigma[i]*z+1)^(1/nu[i]) else y <- mu[i]*exp(sigma[i]*z) 
       if (!is.na(y))  p <- pBCPE(y, mu=mu[i], sigma=sigma[i], nu=nu[i], tau=tau[i])
       else p <- NA
       Y[i] <- y
       P[i] <- p
    }
   Yunder <- ifelse(obj$y < Y, TRUE, FALSE)  
   DF <- data.frame(Y,P,Yunder)
   df <- na.omit(DF)
   sumYu <-  sum(df$Yunder)
   cat("there are ", sumYu, "y values under their minimum turning point \n")
   maxp <- if (sumYu==0) 0 else max(df$P)
   cat("the maximum probability of the lower tail (below minimum turning point) is ", maxp, "\n")
   if (maxp > 0.001) 
     { df1 <- subset(df, df$P>0.001, select=c(Y,P))
      cat("values with high lower tail probability (below minimum turning point) \n")  
      df1
     }  
  }
}
