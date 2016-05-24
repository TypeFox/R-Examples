fun.simu.gld.lm.alt <-
function (n.simu, formula, fit.obj, data, param, fun, maxit = 20000, 
    method = "Nelder-Mead", init = NULL) 
{
if( identical(fun,fun.RMFMKL.ml.m) |
identical(fun,fun.RMFMKL.ml)  |
identical(fun,fun.RPRS.ml.m)  |
identical(fun,fun.RPRS.ml)  ){
optim.fun <- fun.model.lm.optim
}
if( identical(fun,fun.RMFMKL.lm) |
identical(fun,fun.RPRS.lm)  ){
optim.fun <- fun.model.lm.optim.Lmoment
}
    index.var <- 1:(match("L1", names(fit.obj[[3]])) - 1)
    r <- lapply(1:n.simu, function(i, formula, fit.obj, index.var, 
        data, param, fun, maxit, optim.fun, init) {
        print(i)
        GLD.lm.simu(formula, init.coeff = c(fit.obj[[3]][index.var]), 
            init.resid = fit.obj[[3]][(max(index.var) + 1):(max(index.var) + 
                4)], data = data, param = param, fun = fun, method = method, 
            maxit = maxit, optim.fun=optim.fun,fit=fit.obj,init=init)
    }, formula, fit.obj, index.var, data, param, fun, maxit, optim.fun, init)
    r <- do.call("rbind", r)
    return(r)
}
