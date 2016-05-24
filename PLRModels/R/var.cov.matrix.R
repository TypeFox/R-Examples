
###########################################################################
# It fits an optimal ARMA (according to an information criteria) to the data (an only series, that is, a vector).
# It analyses the residual and it returns the variances-covariances matrix of such ARMA model.
# If the ARMA model is not suitable, a message is shown.
###########################################################################

var.cov.matrix <- function(x=1:100, n=4, p.max=3, q.max=3, ic="BIC", 
                           p.arima=NULL, q.arima=NULL, alpha=0.05, num.lb=10) {
# n x n is the dimension of the output matrix

##########
  if (is.null(p.arima) && is.null(q.arima)) {
    p.q <- best.arima(x=x, order.max=c(p.max,0,q.max), include.mean=FALSE, criterio=ic)
    p_opt <- p.q[1,1]
    q_opt <- p.q[1,2]
  }
  
  else if (is.null(p.arima) && !(is.null(q.arima))) {p_opt <- 0; q_opt <- q.arima}
  else if (is.null(q.arima) && !(is.null(p.arima))) {q_opt <- 0; p_opt <- p.arima}
  else {p_opt <- p.arima; q_opt <- q.arima}
##########


# p.q <- best.arima(x=x, order.max=c(p.max,0,q.max), include.mean=FALSE, criterio=ic)[1,]

fitted.model <- arima(x=x, order=c(p_opt,0,q_opt), include.mean=FALSE)

ar.ma <- fitted.model$arma[1:2]

fitdf <- sum(fitted.model$arma[1:2])


pv.lb.t <- c(rep(0,num.lb+1))


for (i in 1:num.lb)
pv.lb.t[i] <- Box.test(x=residuals(fitted.model), lag = fitdf+i, type = "Ljung-Box", fitdf = fitdf)$p.value


pv.lb.t[num.lb + 1] <- t.test(residuals(fitted.model), mu=0)$p.value




if (min(pv.lb.t)<alpha)
	cat("The fitted ARMA model used to estimate V.eps could be not appropriate", "\n")
		

#if ((length(fitted.model$model$phi)+length(fitted.model$model$theta))==0) Var.Cov.x <- diag(1,n)
if (sum(fitted.model$arma[1:2])==0) Var.Cov.x <- diag(1,n)

else
	Var.Cov.x <- toeplitz(ARMAacf(ar=fitted.model$model$phi, ma=fitted.model$model$theta, lag.max=n-1))
  
v.x <- var(x)
v.x <- as.numeric(v.x)

Var.Cov.x <- v.x * Var.Cov.x                 

list(Var.Cov.x=Var.Cov.x, pv.Box.test=pv.lb.t[-(num.lb+1)], pv.t.test=pv.lb.t[(num.lb+1)], ar.ma=ar.ma)

}