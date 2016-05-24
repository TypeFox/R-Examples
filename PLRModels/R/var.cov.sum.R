

###########################################################################
# It fits optimal ARMAs (according to an information criteria) to the data (several series, one in each column of the data matrix).
# It analyses the residual and it returns a vector with the sum of the covariances of such ARMA models.
# If the ARMA model is not suitable, a message is shown.
###########################################################################

var.cov.sum <- function(X=1:100, lag.max=50, p.max=3, q.max=3, ic="BIC", alpha=0.05, num.lb=10) {

X <- as.matrix(X)

L <- ncol(X)

sum.var.cov.X <- 0

for (j in 1:L) {

	x <- X[,j]

	p.q <- best.arima(x=x, order.max=c(p.max,0,q.max), include.mean=FALSE, criterio=ic)[1,]

	fitted.model <- arima(x=x, order=c(p.q[1,1],0,p.q[1,2]), include.mean=FALSE)
  
	ar.ma <- fitted.model$arma[1:2]

	fitdf <- sum(fitted.model$arma[1:2])

  
	pv.lb.t <- c(rep(0,num.lb+1))
  
  
  for (i in 1:num.lb)
	pv.lb.t[i] <- Box.test(x=residuals(fitted.model), lag = fitdf+i, type = "Ljung-Box", fitdf = fitdf)$p.value
  
  pv.lb.t[num.lb + 1] <- t.test(residuals(fitted.model), mu=0)$p.value

  
  
	if (min(pv.lb.t)<alpha)
	  cat("The fitted ARMA model used to estimate Tau.eps could be not appropriate", "\n")
	
	
	 if (sum(fitted.model$arma[1:2])==0) sum.var.cov.X[j] <- 1

		else
			sum.var.cov.X[j] <- 1+2*sum(ARMAacf(ar=fitted.model$model$phi, ma=fitted.model$model$theta, lag.max=lag.max)[-1])
  
			v.x <- var(x)
			v.x <- as.numeric(v.x)
			
			sum.var.cov.X[j] <- v.x * sum.var.cov.X[j] 
			
}

        
list(sum.var.cov.X=sum.var.cov.X, pv.Box.test=pv.lb.t[-(num.lb+1)], pv.t.test=pv.lb.t[(num.lb+1)], ar.ma=ar.ma)

}
