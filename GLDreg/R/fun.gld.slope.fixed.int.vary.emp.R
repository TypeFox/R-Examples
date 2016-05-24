fun.gld.slope.fixed.int.vary.emp <-
function (q, fit, fit.simu, maxit = 20000, method = "Brent") 
{
    
# Use a default lower and upper limit

lubound<-quantile(subset(fit.simu,select=c("(Intercept)")),c(0,1))

lbound<-min(lubound)
ubound<-max(lubound)

    x <- fit$x
    y <- fit$y
    fit <- fit[[3]][-c((length(fit[[3]]) - 3):length(fit[[3]]))]
    k <- apply(fit.simu, 2, function(x, q) quantile(x, q), q)[1]
    r <- optim(k, function(k, x, y, fit, q) {
        resid <- y - data.matrix(x) %*% c(k, fit[-1])
        return((sum(resid <= 0)/length(resid) - q)^2)
    }, x = x, y = y, fit = fit, q = q, control = list(maxit = maxit), 
        method = method,lower=lbound,upper=ubound)
    r.val <- setNames(c(r$par, fit[-1], r$value, r$convergence),c(names(fit),"Objective Value","Convergence"))
    return(list(r, r.val))
}
