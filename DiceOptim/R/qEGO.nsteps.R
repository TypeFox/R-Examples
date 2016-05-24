

qEGO.nsteps <- function (fun, model, npoints, nsteps, lower = rep(0, model@d), upper = rep(1, model@d), crit = "exact", minimization = TRUE,  optimcontrol = NULL, cov.reestim=TRUE, ...) 
{
    history <- rep(NaN,nsteps)
    d <- model@d
    n <- nrow(model@X)
    parinit <- optimcontrol$parinit
    for (i in 1:nsteps) {
	    print(paste("step :",i,"/",nsteps))
	    history[i] <- min(model@y)*minimization + max(model@y)*(!minimization)
	    res <- max_qEI(model=model, npoints = npoints, lower = lower, upper = upper,crit=crit, minimization = minimization, optimcontrol = optimcontrol)
	    model <- update(object = model, newX = res$par, newy=as.matrix(apply(res$par,1, fun, ...), npoints), newX.alreadyExist=FALSE,
                                  cov.reestim = cov.reestim, kmcontrol = list(control = list(trace = FALSE)))
    }
    return(list(par = model@X[(n + 1):(n + nsteps * npoints), 
        , drop = FALSE], value = model@y[(n + 1):(n + nsteps * 
        npoints), , drop = FALSE], npoints = npoints, nsteps = nsteps, 
        lastmodel = model,history=history))
}

