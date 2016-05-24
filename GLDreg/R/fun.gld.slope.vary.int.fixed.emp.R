fun.gld.slope.vary.int.fixed.emp <-
function (q, fit, fit.simu, maxit = 20000, method = "Nelder-Mead") 
{
    x <- fit$x
    y <- fit$y
    k <- apply(fit.simu, 2, function(x, q) quantile(x, q), q)[-1]
    fit <- fit[[3]][-c((length(fit[[3]]) - 3):length(fit[[3]]))]
    r <- optim(k, function(k, x, y, fit, q) {
        resid <- y - data.matrix(x) %*% c(fit[1], k)
        return((sum(resid <= 0)/length(resid) - q)^2)
    }, x = x, y = y, fit = fit, q = q, control = list(maxit = maxit), 
        method = method)
    r.val <- setNames(c(fit[1], r$par, r$value, r$convergence),c(names(fit),"Objective Value","Convergence"))
    return(list(r, r.val))
}
