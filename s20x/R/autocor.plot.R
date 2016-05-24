autocor.plot<- function(fit){
    current.res <- fit$residuals[-1]
    T <- length(fit$residuals)
    lagged.res <- fit$residuals[-T]
    plot(lagged.res, current.res, main = "Current vs Lagged residuals")
ab<-lm(current.res~lagged.res)$coeff
##abline(ab[1],ab[2],lty=2)
abline(v=0,lty=3)
abline(h=0,lty=3)
}
