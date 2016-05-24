
plrm.ci <- function(data=data, seed=123, CI="AD", B=1000, N=50, a=NULL,
                    b1=NULL, b2=NULL,
                    estimator="NW", kernel="quadratic",
                    p.arima=NULL, q.arima=NULL, p.max=3, q.max=3, 
                    alpha=0.05, alpha2=0.05, num.lb=10, ic="BIC", Var.Cov.eps=NULL) 
{
  
  if (!is.matrix(data))  stop("data must be a matrix")
  if (ncol(data) < 3)  stop("data must have at least 3 columns")
  
  if (!is.numeric(seed)) stop ("seed must be numeric") 
  set.seed(seed)
  
  if ((CI != "AD") & (CI != "B")  & (CI != "all"))  stop("IC=AD or IC=B or IC=all is required")
  
  if (is.null(B))   stop ("B must not be NULL") 
  if (length(B) !=1)  stop ("B must be an only value")
  if (!is.numeric(B))   stop ("B must be numeric") 
  if (B<1)  stop ("B must be a positive value") 
  
  if (is.null(N))   stop ("N must not be NULL") 
  if (length(N) !=1)  stop ("N must be an only value")
  if (!is.numeric(N))   stop ("N must be numeric") 
  if (N<0)  stop ("N must be a positive value") 
  
  if ((!is.null(b1)) && (length(b1) !=1))  stop ("b1 must be an only value")
  if ((!is.null(b1)) && (!is.numeric(b1)))   stop ("b1 must be numeric") 
  if ((!is.null(b1)) && (b1<0))  stop ("b1 must be a positive value") 
  
  if ((!is.null(b2)) && (length(b2) !=1))  stop ("b2 must be an only value")
  if ((!is.null(b2)) && (!is.numeric(b2)))   stop ("b2 must be numeric") 
  if ((!is.null(b2)) && (b2<0))  stop ("b2 must be a positive value") 
  
  if ((estimator != "NW") & (estimator != "LLP"))  stop("estimator=NW or estimator=LLP is required")
  
  if ((kernel != "quadratic") & (kernel != "Epanechnikov") & (kernel != "triweight") & (kernel != "gaussian") & (kernel != "uniform"))  stop("kernel must be one of the following: quadratic, Epanechnikov, triweight, gaussian or uniform")
  
  if ((!is.null(p.arima)) && (length(p.arima) !=1))  stop ("p.arima must be an only value")
  if ((!is.null(p.arima)) && (!is.numeric(p.arima)))   stop ("p.arima must be numeric") 
  if ((!is.null(p.arima)) && (p.arima<0))  stop ("p.arima must be a positive value") 
  
  if ((!is.null(q.arima)) && (length(q.arima) !=1))  stop ("q.arima must be an only value")
  if ((!is.null(q.arima)) && (!is.numeric(q.arima)))   stop ("q.arima must be numeric") 
  if ((!is.null(q.arima)) && (q.arima<0))  stop ("q.arima must be a positive value") 
  
  if (is.null(p.max))   stop ("p.max must not be NULL") 
  if (length(p.max) !=1)  stop ("p.max must be an only value")
  if (!is.numeric(p.max))   stop ("p.max must be numeric") 
  if (p.max<0)  stop ("p.max must be a positive value") 
  
  if (is.null(q.max))   stop ("q.max must not be NULL") 
  if (length(q.max) !=1)  stop ("q.max must be an only value")
  if (!is.numeric(q.max))   stop ("q.max must be numeric") 
  if (q.max<0)  stop ("q.max must be a positive value") 
  
  if (is.null(alpha))   stop ("alpha must not be NULL") 
  if (length(alpha) !=1)  stop ("alpha must be an only value")
  if (!is.numeric(alpha))   stop ("alpha must be numeric") 
  if ( (alpha<0) | (alpha>1) )  stop ("alpha must be between 0 and 1") 
  
  if (is.null(alpha2))   stop ("alpha2 must not be NULL") 
  if (length(alpha2) !=1)  stop ("alpha2 must be an only value")
  if (!is.numeric(alpha2))   stop ("alpha2 must be numeric") 
  if ( (alpha2<0) | (alpha2>1) )  stop ("alpha2 must be between 0 and 1") 
  
  if (is.null(num.lb))   stop ("num.lb must not be NULL") 
  if (length(num.lb) !=1)  stop ("num.lb must be an only value")
  if (!is.numeric(num.lb))   stop ("num.lb must be numeric") 
  if (num.lb<=0)  stop ("num.lb must be a positive value") 
  
  if ( (ic != "BIC") & (ic != "AIC") & (ic != "AICC") )  stop("ic=BIC or ic=AIC or ic=AICC is required")
  
  if ( (!is.null(Var.Cov.eps)) && (sum(is.na(Var.Cov.eps))  != 0) ) stop("Var.Cov.eps must have numeric values")
  if ( (!is.null(Var.Cov.eps)) && (!is.matrix(Var.Cov.eps)) )  stop("Var.Cov.eps must be a matrix")
  if ( (!is.null(Var.Cov.eps)) && ( (ncol(Var.Cov.eps) != nrow(data)) | (nrow(Var.Cov.eps) != nrow(data)) ) ) stop("Var.Cov.eps must have dimension n x n")
  if ( (!is.null(Var.Cov.eps)) && (any(t(Var.Cov.eps) != Var.Cov.eps)  ) ) stop("Var.Cov.eps must be symmetric")
  
  # ###########################################################
  # ###########################################################
  
  n <- nrow(data)  
  p <- ncol(data)-2
  
  Y <- data[, 1]
  X <- data[, 2:(p+1)]
  t <- data[, p+2]
  
  if (!is.matrix(X))  X <- as.matrix(X) 
  
  if (is.null(a)) {
    if (p==1) a <- c(1)
    else a <- c(1,rep(0,p-1))
  }
  else if (!is.vector(a))   stop ("a must be a vector") 
  else if (!is.numeric(a))   stop ("a must be numeric") 
  else if (p!=length(a)) stop("vector a must have length: ncol(data)-1")
  
  if (is.null(b1)) 
    b1 <- plrm.cv(data=data, estimator=estimator, kernel=kernel)$bh.opt[2,1]
  
  
  PLRM <- plrm.est(data, b=b1, estimator=estimator, kernel=kernel)
  beta.est <- PLRM$beta
  eps <- PLRM$residuals
  
  
  if (is.null(p.arima) && is.null(q.arima)) {
    p.q <- best.arima(x=eps, order.max=c(p.max,0,q.max), include.mean=FALSE)
    p_opt <- p.q[1,1]
    q_opt <- p.q[1,2]
  }
  
  
  else if (is.null(p.arima) && !(is.null(q.arima))) {p_opt <- 0; q_opt <- q.arima}
  else if (is.null(q.arima) && !(is.null(p.arima))) {q_opt <- 0; p_opt <- p.arima}
  else {p_opt <- p.arima; q_opt <- q.arima}
  
  
  eps.arima <- arima(x=eps, order=c(p_opt,0,q_opt), include.mean=FALSE)
  
  
  if (CI=="B") {
    
    b.pv.t <- c(rep(NA,num.lb+1))
    fitdf <- sum(eps.arima$arma[1:2])
    
    for (i in 1:num.lb)
      b.pv.t[i] <- Box.test(x=residuals(eps.arima), lag = fitdf+i, type = "Ljung-Box", fitdf = fitdf)$p.value
    b.pv.t[i+1] <- t.test(residuals(eps.arima), mu=0)$p.value
    
    if (min(b.pv.t)<alpha2)
      cat("The fitted ARMA model could be not appropriate", "\n")
    
  }
  
  
  if ((CI=="AD") | (CI=="all")) {
    
    if (is.null(Var.Cov.eps)) {
      Var.Cov.eps <- matrix(NA, n, n)
      Var.Cov.mat <- var.cov.matrix(x=eps, n=n, p.max=p.max, q.max=q.max, ic=ic, p.arima=p_opt, q.arima=q_opt, alpha=alpha, num.lb=num.lb)
      
      Var.Cov.eps <- Var.Cov.mat[[1]]
      ad.pv.Box.test <- Var.Cov.mat[[2]]
      ad.pv.t.test <- Var.Cov.mat[[3]]
    }
    
    XX <- matrix(0,n,p)
    WX <- matrix(0,n,1)
    
    
    for (j in 1:p) {
      data1 <- cbind(data[,j+1],t)
      WX <- np.est(data=data1,h.seq=b1,newt=t,estimator=estimator,kernel=kernel)
      XX[,j] <- data[,j+1]-WX
    }
    
    
    X.X <- t(XX)%*%XX
    X.X.1 <- solve(X.X)
    
    A <- n*X.X.1%*%t(XX)%*%Var.Cov.eps%*%XX%*%X.X.1
    
    Dt <- sqrt((t(a)%*%A%*%a)/n)
    
    z.quantile <- qnorm(1-alpha/2)
    
    
    ci_inf_ad <- t(a)%*%beta.est - z.quantile * Dt
    ci_sup_ad <- t(a)%*%beta.est + z.quantile * Dt
    
  }
  
  if ((CI=="B") | (CI=="all")) {
    if (p_opt==0) ar.coef <- 0
    else ar.coef <- as.numeric(eps.arima$coef[1:p_opt])
    
    if (q_opt==0) ma.coef <- 0
    else ma.coef <- as.numeric(eps.arima$coef[(p_opt+1):(p_opt+q_opt)])
    
    
    white.noise <- eps.arima$residuals
    white.noise <- white.noise[(p_opt+1):n]
    
    white.noise <- white.noise-mean(white.noise)
    
    beta.est.boots <- matrix(NA,B,p)
    for (m in 1:B) {
      
      
      white.noise.sample <- sample(white.noise, size=n+N, replace=TRUE)
      
      fi <- c(rep(NA,N+1))     
      fi[1] <- 1
      fi[-1] <- ARMAtoMA(ar=ar.coef, ma=ma.coef, lag.max=N)
      fi <- rev(fi)
      
      eps.boots <- c(rep(NA,n))
      for (i in 1:n) {
        eps.ij.boots <- white.noise.sample[i:(i+N)]
        
        eps.boots[i] <-  fi%*%eps.ij.boots
        
      }
      
      
      Y.boots <- PLRM$fitted.values + eps.boots
      
      data.boots <- cbind(Y.boots, X, t)
      
      if ((m==1) & (is.null(b2)))
        b2 <- plrm.cv(data=data.boots, estimator=estimator, kernel=kernel)$bh.opt[2,1]
      
      beta.est.boots[m,] <- plrm.est(data.boots, b=b2, estimator=estimator, kernel=kernel)$beta        
      
    }
    
    beta.est <- as.vector(beta.est)
    dif.beta <- t(beta.est.boots)-beta.est
    z <- sqrt(n)*t(a)%*%dif.beta
    
    z.quantile1 <- quantile(z,1-alpha/2)
    z.quantile2 <- quantile(z,alpha/2)
    
    ci_inf_b <- a%*%beta.est - (1/sqrt(n))*z.quantile1
    ci_sup_b <- a%*%beta.est - (1/sqrt(n))*z.quantile2
    
  }
  
  if (CI=="B") list(Bootstrap=data.frame(ci_inf=ci_inf_b, ci_sup=ci_sup_b, p_opt=p_opt, q_opt=q_opt, b1=b1, b2=b2), pv.Box.test=b.pv.t[-num.lb+1], pv.t.test=b.pv.t[num.lb+1])
  else if (CI=="AD") list(AD=data.frame(ci_inf=ci_inf_ad, ci_sup=ci_sup_ad, p_opt=p_opt, q_opt=q_opt, b1=b1), pv.Box.test=ad.pv.Box.test, pv.t.test=ad.pv.t.test)
  else if (CI=="all") list(Bootstrap=data.frame(ci_inf=ci_inf_b, ci_sup=ci_sup_b, p_opt=p_opt, q_opt=q_opt, b1=b1, b2=b2),AD=data.frame(ci_inf=ci_inf_ad, ci_sup=ci_sup_ad, p_opt=p_opt, q_opt=q_opt, b1=b1), pv.Box.test=ad.pv.Box.test, pv.t.test=ad.pv.t.test)
  
}

