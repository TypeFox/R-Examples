fun.gld.all.vary.emp <-
function (q, fit, fit.simu, maxit = 20000, method = "Nelder-Mead") 
{
    x <- fit$x
    y <- fit$y
    k <- apply(fit.simu, 2, function(x, q) quantile(x, q), q)
    r <- optim(k, function(k, x, y, q) {
        resid <- y - x %*% c(k)
        return((sum(resid <= 0)/length(resid) - q)^2)
    }, x = x, y = y, q = q, control = list(maxit = maxit), method = method)
    r.val <- setNames(c(r$par, r$value, r$convergence),c(names(fit$"Estimated"[-c((length(fit$"Estimated")-3):length(fit$"Estimated"))]),"Objective Value","Convergence"))
    return(list(r, r.val))
}
