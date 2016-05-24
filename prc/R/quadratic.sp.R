# When fix.sigma.sq is TRUE, 
quadratic.eiv.sp=function(xvar, yvar, grid.density=200, init=NULL, reltol=1e-3, opt.method=c("optim"), max.iter=50, fix.sigma.sq=FALSE, verbose=FALSE) {    
    
    # removing this line requires adding package namespace to Matrix and mosek function calls
    have.Rmosek=requireNamespace("Rmosek")
    if (!have.Rmosek) {
            print("Rmosek does not load successfully")
            return (NULL)
    }
    
    stopifnot(length(xvar)==length(yvar))
    n=length(xvar)    
    xvar=as.vector(xvar); yvar=as.vector(yvar) # strip off potential attributes; otherwise it may cause trouble with gnls
    dat=data.frame(readout.x=unname(xvar), readout.y=unname(yvar))    
    
    res=list(xvar=xvar, yvar=yvar)
    class(res)=c("quad",class(res))
    
    #### run quad for a small number of iterations to initialize
    if (is.null(init)) {
        if (verbose) cat("run quadratic.eiv to get initial theta\n")
        fit.init=quadratic.eiv(xvar, yvar, verbose=verbose)
        theta=c(coef(fit.init), sigma.sq=fit.init$sigma.sq)
    } else {
        theta=init
        if (length(theta)!=4) stop("init need to include sigma.sq")
    }
    
    #### Alternate between two steps
    
    if (verbose) cat("\nInitial theta:", theta, "\n")
    iterations=0
    
    while(TRUE) {
    
        iterations=iterations+1
        
        # Step 2: get nonparametric estimate of distribution of u
        
        # define the grid and compute A
        u=seq(min(xvar), max(xvar), length=2+grid.density)
        u=u[-c(1,grid.density)] # length of u is K
        A=compute.A.quad (theta["a"], theta["b"], theta["c"], theta["sigma.sq"], u, xvar, yvar)         
        
        # run mosek
        sco1 <- list(sense = "min")
        sco1$c <- rep(0,n)
        sco1$A <- Matrix::Matrix (t(A), sparse = TRUE )
        sco1$bc <- rbind(blc = rep(-Inf,length(u)), buc = rep(n,length(u)))
        sco1$bx <- rbind(blx = rep(0,n), bux = rep(Inf,n))
        opro <- matrix(list(), nrow=5, ncol=n, dimnames=list(c("type","j","f","g","h"), NULL))
        for (i in 1:n) opro[,i] <- list("LOG", i, -1.0, 1.0, 0.0)
        sco1$scopt <- list(opro=opro)
        fit.mosek <- Rmosek::mosek(sco1, opts=list(verbose=1))
        if (fit.mosek$sol$itr$solsta!="OPTIMAL") warning("fit.mosek$sol$itr$solsta!=OPTIMAL")
        if (fit.mosek$sol$itr$prosta!="PRIMAL_AND_DUAL_FEASIBLE") warning("fit.mosek$sol$itr$prosta!=PRIMAL_AND_DUAL_FEASIBLE")
    
        if (verbose>=3) {
            plot(fit.mosek$sol$itr$xc, sco1$A %*% fit.mosek$sol$itr$xx); abline(0,1) # sanity check
            print(table(fit.mosek$sol$itr$skc))        
            #print(fit.mosek)
#            print(table(fit.mosek$sol$itr$skx))
#            print(summary(fit.mosek$sol$itr$xx))
#            print(summary(fit.mosek$sol$itr$xc))
#            print(cbind(fit.mosek$sol$itr$skx, fit.mosek$sol$itr$xx, fit.mosek$sol$itr$slx, fit.mosek$sol$itr$sux))
            print(cbind(fit.mosek$sol$itr$skc, fit.mosek$sol$itr$xc, fit.mosek$sol$itr$slc, fit.mosek$sol$itr$suc))
        }
    
        # recover primal soln
        #support.set = which(fit.mosek$sol$itr$skc=="UL") # UL seems to correspond to a much higher threshold, something like 0.1, also when this is used, not getting reproducible results
        support.set = which(fit.mosek$sol$itr$suc>0.001) 
        if (length(support.set)==0) { print("no support set"); break}
        p.fit=lm.fit(x = A[,support.set], y = 1/fit.mosek$sol$itr$xx)
        if (verbose) plot(p.fit$fitted.values, 1/fit.mosek$sol$itr$xx) # making sure we are doing a good job at recovering primal soln
        p=coef(p.fit)
        support=u[support.set]
        if (any(is.na(p))) { 
            p[is.na(p)]=-1
            #stop("some p are NA")
        }
        mix.lik = sum(-log(fit.mosek$sol$itr$xx))         
        if (verbose) {
            if (verbose>=3) {
                mix.lik.1 = sum(log(compute.A.quad (theta["a"], theta["b"], theta["c"], theta["sigma.sq"], support, xvar, yvar) %*%p))
                print(c(mix.lik, mix.lik.1))
                print(rbind(support, p))
            }
            plot(0,0, type="n", ylim=c(0,max(p)), xlim=range(support), xlab="r", ylab="p", main="Iteration "%+%iterations)
            points(support, p, pch=19)
            for (i in 1:length(p)) lines(rep(support[i],2), c(0,p[i]))            
        }
        # sometimes, some p are negative
        while (any(p<0)) {
            if (verbose) print("some p smaller than 0")
            support.set=support.set[p>0]
            p.fit=lm.fit(x = A[,support.set], y = 1/fit.mosek$sol$itr$xx)
            p=coef(p.fit)
            support=u[support.set]
            if (verbose>=3) print(rbind(support, p))
            if (verbose) {
                plot(0,0, type="n", ylim=c(0,max(p)), xlim=range(support), xlab="r", ylab="p", main="Iteration "%+%iterations)
                points(support, p, pch=19)
                for (i in 1:length(p)) lines(rep(support[i],2), c(0,p[i]))            
            }
        }
    
        # Step 3: estimate theta and sigma.sq
        # don't think optim BFGS results depend on seed
        if (!fix.sigma.sq) {
            optim.out = optim(
                par=theta,
                fn = function(theta,...) {
                    a=theta["a"]; b=theta["b"]; c=theta["c"]; sigma.sq=theta["sigma.sq"]
                    A=compute.A.quad (a,b,c, sigma.sq, support, xvar, yvar) 
                    -mean(log(A%*%p))
                }, 
                gr=NULL,
                p, support, xvar, yvar, 
                method="BFGS", control = list(trace=0), hessian = F
            )
            new.theta = optim.out$par
            if (verbose>=4) print(optim.out)
            mix.lik.2 = sum(log(compute.A.quad (new.theta["a"],new.theta["b"],new.theta["c"],new.theta["sigma.sq"], support, xvar, yvar) %*%p))
        } else {
            optim.out = optim(
                par=theta,
                fn = function(theta,...) {
                    a=theta["a"]; b=theta["b"]; c=theta["c"]; #sigma.sq=theta["sigma.sq"]
                    A=compute.A.quad (a,b,c, sigma.sq=init["sigma.sq"], support, xvar, yvar) 
                    -mean(log(A%*%p))
                }, 
                gr=NULL,
                p, support, xvar, yvar, 
                method="BFGS", control = list(trace=0), hessian = F
            )
            new.theta = optim.out$par
            if (verbose>=4) print(optim.out)
            mix.lik.2 = sum(log(compute.A.quad (new.theta["a"],new.theta["b"],new.theta["c"],init["sigma.sq"], support, xvar, yvar) %*%p))
        }
        
            
        if (verbose) {
            cat("Iter "%+%iterations%+%".", "UL", length(support), "mix.lik", mix.lik, "mix.lik.2", mix.lik.2, "theta:", new.theta)#, "support:", mean(support))
            cat("\n")            
        }
        if (max(abs(1 - new.theta/theta)) < reltol) {
            if (verbose) cat("converged\n")
            theta = new.theta # update theta 
            break;
        } else if (iterations>=max.iter) {
            if (verbose) cat("max iter reached\n")
            theta = new.theta # update theta 
            break;
        } else {
            theta = new.theta # update theta 
        }        
            
    } # end while loop
        
    res$support=support
    res$p=p
    res$coefficients=theta[1:3]
    res$sigma.sq = ifelse (!fix.sigma.sq, theta["sigma.sq"], init["sigma.sq"])
    res$A=with(res, compute.A.quad (coefficients["a"],coefficients["b"],coefficients["c"], sigma.sq, support, xvar, yvar) )
    
    res
}

mixlik.quad=function(fit) {    
    with(fit, sum(log(compute.A.quad (coefficients["a"],coefficients["b"],coefficients["c"],sigma.sq, support, xvar, yvar) %*%p)))
    
}

# populate A
compute.A.quad = function(a,b,c, sigma.sq, support, xvar, yvar) {
    .Call("compute_A_quad", a,b,c, sigma.sq, support, xvar, yvar)
}
