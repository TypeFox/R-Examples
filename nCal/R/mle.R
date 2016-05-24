# fit.4pl=FALSE; var.model="power"; robust="mean"; method="drm"; max.iter=20; reltol=1e-3; verbose=TRUE; log.both.sides=F
mle.dr=function (formula, data, start, max.iter=50, reltol=1e-3, fit.4pl=FALSE, 
    verbose=FALSE) {


    outcome.var=model.frame(formula, data)[[1]]
    data[["outcome.var"]]=outcome.var # this is necessary, b/c drc uses model.weights to extract weight from a model.frame
    predictor.coln=all.vars(formula)[2]    
    
    # iterate between drm fit given a fixed weights and estimating the variance component parameters
    iterations=0
    theta=start[c("c","d","g","h","f")]
    g.hat=start["g.hat"]
    sigmasq.hat=start["sigmasq.hat"]
    devs=c()

    f1=function(x,c,d,g,h,f) (c + (d - c)/(1 + exp(-log(f) - h/(d - c) * (1 + 1/f)^(1 + f) * (log(x) - g)))^f)
    fit=list()
    while (TRUE) {
        iterations=iterations+1
        
        # estimate curve fit via optim
        optim.out= optim(
            par=theta, 
            function(theta.fn,...) {
                fitte=f1(data[[predictor.coln]], theta.fn["c"],theta.fn["d"],theta.fn["g"],theta.fn["h"],ifelse(fit.4pl,1,theta.fn["f"]))
                resi = outcome.var-fitte
                sum( g.hat*log(fitte) + resi^2/fitte^g.hat/sigmasq.hat )
            }, 
            gr=NULL,
            method="BFGS", control = list(trace=0), hessian = F
        )             
        new.theta=optim.out$par
            
        # estimate variance component, g.hat and sigmasq.hat, via uniroot
        fitte=f1(data[[predictor.coln]], new.theta["c"],new.theta["d"],new.theta["g"],new.theta["h"],ifelse(fit.4pl,1,new.theta["f"]))
        resi = outcome.var-fitte        
        new.g.hat = try(uniroot(f=function(g,...) sum( 
            log(fitte) * ( 1 - resi^2/{fitte^g * mean(resi^2/fitte^g)} )
        ), interval=c(-3,3), fitte)$root, silent=TRUE) 
        # try one more time if an error happens, this time widen the interval
        if (inherits(new.g.hat, "try-error")) {
            new.g.hat = try(uniroot(f=function(g,...) 
            sum( 
                log(fitte) * ( 1 - resi^2/{fitte^g * mean(resi^2/fitte^g)} ) # mean(resi^2/fitte^g) is sigmasq.hat
            ), interval=c(-10,10), fitte)$root) 
            if (inherits(new.g.hat, "try-error")) {
                cat("fail to find gamma, return null\n")
                return (NULL)
            }
        }
        new.sigmasq.hat = mean(resi^2/fitte^new.g.hat)
        
        # compute deviance
        new.dev=sum( new.g.hat*log(fitte) + log(new.sigmasq.hat) + resi^2/(fitte^new.g.hat*new.sigmasq.hat) )  
        devs = c(devs, new.dev)
        
        if (verbose) cat("Iter "%+%iterations%+%".", "theta:", new.theta, "g.hat:", new.g.hat, " deviance:", devs[iterations], "\n")
        
        # stopping rules
        if (max(abs(1 - new.theta/theta)) < reltol) {
            if (verbose) cat("converged\n")
            fit$converged=TRUE
            break;
        } else if (iterations>=max.iter) {
            if (verbose) cat("max iter reached\n")
            fit$converged=FALSE
            break;
        } else if (devs[iterations]>ifelse(iterations>1,devs[iterations-1],Inf)) {
            if (verbose) cat("dev increasing, stop and revert to previous fit\n")
            fit$converged=FALSE
            break;
        } else {
            # update theta, save fit to old.fit, and continue
            theta = new.theta 
            g.hat=new.g.hat
            sigmasq.hat=new.sigmasq.hat
            dev = new.dev
        }        
        
    } # end loop
    
    fit$coefficients=theta
    fit$var.power = g.hat
    fit$sigmasq.hat = sigmasq.hat
    fit$deviance=dev
    fit$iter=iterations
    
    fit
    
}
