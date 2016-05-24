#' Model selection at a known breakpoint
#' 
#' Returns all parameter estimates and log-likelihoods for all possible models at a selected breakpoint.  These are: \itemize{ \item{M0} - all parameters equal \item{M1} - \eqn{\mu_1 != \mu_2} \item{M2} - \eqn{\sigma_1 != \sigma_2} \item{M3} - \eqn{\tau_1 != \tau_2} \item{M4} - \eqn{\mu_1 != \mu_2} AND \eqn{\sigma_1 != \sigma_2} \item{M5} - \eqn{\mu_1 != \mu_2} AND \eqn{\tau_1 != \tau_2} \item{M6} - \eqn{\sigma_1 != \sigma_2} AND \eqn{\tau_1 != \tau_2} \item{M7} - all parameters unequal}

#' @param x  vector of time series values.
#' @param t  vector of times of measurements associated with x.
#' @param breakpoint breakpoint to test (in terms of the INDEX within "t" and "x", not actual time value).
#' @param K sensitivity parameter.  Standard definition of BIC, K = 2.  Smaller values of K mean less sensitive selection, i.e. higher likelihood of selecting null (or simpler) models. 
#' @param tau whether or not to estimate time scale \eqn{\tau} (preferred) or autocorrelation \eqn{\rho}.

#' @return Returns a names matrix with 8 rows (one for each model) and columns: \code{Model}, \code{LL, bic, mu1, s1, rho1, mu2, s2, rho2}.  Fairly self-explanatory.  Note that the \code{rho} columns include the \code{tau} values, if \code{tau} is TRUE (as it usually should be).
#' 
#' @seealso Used directly within \code{\link{WindowSweep}}.  Relies heavily on \code{\link{GetRho}}.
#' @author Eliezer Gurarie
 
GetModels <-
  function(x,t,breakpoint,K=2, tau = TRUE)
  {
    
    GetLL <-
      function(x,t,mu,s,rho,tau)
      {
        dt <- diff(t)
        n<-length(x)
        x.plus <- x[-1]
        x.minus <- x[-length(x)]
        
        if(tau) RhoPower <- exp(-dt/rho)
        else RhoPower <- rho^dt
        
        Likelihood <- dnorm(x.plus,mean=mu+(RhoPower)*(x.minus-mu),sd=s*sqrt(1-RhoPower^2))
        
        LL <- sum(log(Likelihood))
        return(LL)
      }
    
    
    M0 <-
      function(x,t,breakpoint,K=2, ...)
        # null model: all mus, s's, rhos the same
      {
        rhoLL <- GetRho(x,t, ...)
        LL <- rhoLL[2]
        
        bic <- -K*LL + 3*log(length(x))
        
        rho1<-rhoLL[1]
        rho2<-rho1
        
        mu1<-mean(x)
        mu2<-mu1
        
        s1<-sd(x)
        s2<-s1
        
        return(c(LL,bic,mu1,s1,rho1,mu2,s2,rho2))
      }
    
    
    M1 <-
      function(x,t,breakpoint,K=2,...)
        # mus different, all else the same
      {
        
        x1<-x[1:breakpoint]
        x2<-x[(breakpoint+1):length(x)]
        t1<-t[1:breakpoint]
        t2<-t[(breakpoint+1):length(x)]
        
        mu1<-mean(x1)
        mu2<-mean(x2)
        
        xprime <- c(x1-mu1,x2-mu2)
        s1<-sd(xprime)
        s2<-s1
        
        rho1<-as.numeric(GetRho(xprime,t,...)[1])
        rho2<-rho1
        
        LL1<-GetLL(x1,t1,mu1,s1,rho1,...)
        LL2<-GetLL(x2,t2,mu2,s2,rho2,...)
        LL<-LL1+LL2
        bic <- -K*LL + 5*log(length(x))
        
        return(c(LL,bic,mu1,s1,rho1,mu2,s2,rho2))
      }
    
    
    M2 <-
      function(x,t,breakpoint,K=2,...)
        # sds different, all else same
      {
        
        x1<-x[1:breakpoint]
        x2<-x[(breakpoint+1):length(x)]
        
        t1<-t[1:breakpoint]
        t2<-t[(breakpoint+1):length(x)]
        
        mu1<-mean(x)
        mu2<-mu1
        
        s1<-sd(x1)
        s2<-sd(x2)
        
        xprime <- c( (x1-mu1)/s1 , (x2-mu2)/s2 )
        rho1 <- as.numeric(GetRho(xprime,t, ...)[1])
        rho2 <- rho1
        
        LL1<-GetLL(x1,t1,mu1,s1,rho1,...)
        LL2<-GetLL(x2,t2,mu2,s2,rho2,...)
        LL<-LL1+LL2
        bic <- -K*LL + 5*log(length(x))
        
        return(c(LL,bic,mu1,s1,rho1,mu2,s2,rho2))
      }
    
    
    M3 <-
      function(x,t,breakpoint,K=2,...)
        # rhos different, all else same
      {
        x1<-x[1:breakpoint]
        x2<-x[(breakpoint+1):length(x)]
        t1<-t[1:breakpoint]
        t2<-t[(breakpoint+1):length(x)]
        
        mu1 <- mean(x)
        mu2 <- mu1
        
        s1 <- sd(x)
        s2 <- s1
        
        rho1<-as.numeric(GetRho(x1,t1,...)[1])
        rho2<-as.numeric(GetRho(x2,t2,...)[1])
        
        LL1<-GetLL(x1,t1,mu1,s1,rho1,...)
        LL2<-GetLL(x2,t2,mu2,s2,rho2,...)
        LL<-LL1+LL2
        bic <- -K*LL + 5*log(length(x))
        
        return(c(LL,bic,mu1,s1,rho1,mu2,s2,rho2))
      }
    
    
    M4 <-
      function(x,t,breakpoint,K=2,...)
        # mu and sigma different, rho same
      {
        x1<-x[1:breakpoint]
        x2<-x[(breakpoint+1):length(x)]
        t1<-t[1:breakpoint]
        t2<-t[(breakpoint+1):length(x)]
        
        mu1<-mean(x1)
        mu2<-mean(x2)
        s1<-sd(x1)
        s2<-sd(x2)
        
        xprime <- c( (x1-mu1)/s1 , (x2-mu2)/s2 )
        rho1 <- as.numeric(GetRho(xprime,t,...)[1])
        rho2 <- rho1
        
        LL1<-GetLL(x1,t1,mu1,s1,rho1,...)
        LL2<-GetLL(x2,t2,mu2,s2,rho2,...)
        LL<-LL1+LL2
        bic <- -K*LL + 6*log(length(x))
        
        return(c(LL,bic,mu1,s1,rho1,mu2,s2,rho2))
      }
    
    
    M5 <-
      function(x,t,breakpoint,K=2,...)
        # mu and rho different, sigma same
      {
        x1<-x[1:breakpoint]
        x2<-x[(breakpoint+1):length(x)]
        t1<-t[1:breakpoint]
        t2<-t[(breakpoint+1):length(x)]
        
        mu1<-mean(x1)
        mu2<-mean(x2)
        
        xprime <- c(x1-mu1, x2-mu2)
        s1<-sd(xprime)
        s2<-s1
        
        rho1<-as.numeric(GetRho(x1,t1,...)[1])
        rho2<-as.numeric(GetRho(x2,t2,...)[1])
        
        LL1<-GetLL(x1,t1,mu1,s1,rho1,...)
        LL2<-GetLL(x2,t2,mu2,s2,rho2,...)
        LL<-LL1+LL2
        bic <- -K*LL+ 6*log(length(x))
        
        return(c(LL,bic,mu1,s1,rho1,mu2,s2,rho2))
      }
    
    M6 <-
      function(x,t,breakpoint,K=2,...)
        # sigma and rho different, mu same
      {
        x1<-x[1:breakpoint]
        x2<-x[(breakpoint+1):length(x)]
        t1<-t[1:breakpoint]
        t2<-t[(breakpoint+1):length(x)]
        
        mu1<-mean(x)
        mu2<-mean(x)
        s1<-sd(x1)
        s2<-sd(x2)
        
        x1prime <- (x1-mu1)/s1 
        x2prime <- (x2-mu2)/s2
        
        rho1 <- as.numeric(GetRho(x1prime,t1,...)[1])
        rho2 <- as.numeric(GetRho(x2prime,t2,...)[1])
        
        LL1<-GetLL(x1,t1,mu1,s1,rho1,...)
        LL2<-GetLL(x2,t2,mu2,s2,rho2,...)
        
        LL<-LL1+LL2
        
        bic <- -K*LL + 6*log(length(x))
        
        return(c(LL,bic,mu1,s1,rho1,mu2,s2,rho2))
      }
    
    M7 <-
      function(x,t,breakpoint,K=2,...)
        # most "alternative" model: all mus, s's, rhos different
      {
        rhoLL1 <- GetRho(x[1:breakpoint],t[1:breakpoint],...)
        rhoLL2 <- GetRho(x[(breakpoint+1):length(x)],t[(breakpoint+1):length(x)],...)
        
        LL1 <- rhoLL1[2]
        LL2 <- rhoLL2[2]
        
        x1<-x[1:breakpoint]
        x2<-x[(breakpoint+1):length(x)]
        t1<-t[1:breakpoint]
        t2<-t[(breakpoint+1):length(x)]
        
        mu1<-mean(x1)
        mu2<-mean(x2)
        s1<-sd(x1)
        s2<-sd(x2)
        rho1 <- rhoLL1[1]
        rho2 <- rhoLL2[1]
        
        LL <- LL1+LL2
        bic <- -K*LL + 7*log(length(x))
        return(c(LL,bic,mu1,s1,rho1,mu2,s2,rho2))
      }
    
    
    
    r <- matrix(NA, nrow=8, ncol = 8)
    colnames(r) <- c("LL", "bic", "mu1", "s1", "rho1", "mu2", "s2", "rho2")
    
    for(i in 1:8)
    {
      f <- get(paste("M",i-1,sep=""))
      r[i,] <- f(x,t,breakpoint,K, tau)
    }
    
    return(cbind(Model = 0:7, r))
  }
