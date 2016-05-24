logdbf=c("logd","b","f")

# dil.r is dilution ratio, previously named k


# verbose=1: only one line per iteration
# verbose=2: print after step 2, step 3a, step 3b
# verbose=4: print optim output
#reltol=1e-3; max.iter=20; verbose=T; model="sp"; grid.density=200; method.init="TLS"; init=NULL
prcsp=function(xvar, dil.x, yvar, dil.y, 
    model=c("sp","struct"), #sp: semiparametric, struct: structural
    stop.when.dropping=FALSE, grid.density=200, 
    init=NULL, method.init=c("TLS","naive"), reltol=1e-3, max.iter=20, 
    try.additiona.support.sets=FALSE,  # sp only 
    keep.history=FALSE, verbose=FALSE) 
{    
    
    method.init=match.arg(method.init)
    model=match.arg(model)
    if(verbose) myprint(model, method.init)
    
    # removing this line requires adding package namespace to Matrix and mosek function calls
    if(model=="sp") {
        have.Rmosek=requireNamespace("Rmosek")
        if (!have.Rmosek) {
            print("Rmosek does not load successfully")
            return (NULL)
        }
    }
                
    dil.r=dil.x/dil.y
    stopifnot(length(xvar)==length(yvar))
    n=length(xvar)    
    xvar=as.vector(xvar); yvar=as.vector(yvar) # strip off potential attributes; otherwise it may cause trouble with gnls
    dat=data.frame(readout.x=unname(xvar), readout.y=unname(yvar))    
    
    res=list(dilution.ratio=dil.r, dilution.x=dil.x, dilution.y=dil.y, xvar=xvar, yvar=yvar)
    class(res)=c("prc",class(res))
    
    #########################################################
    # Step 1: run prc to the end to initialize
    
    if (is.null(init)) {
        if (verbose) cat("\n========== run prc to initialize ================\n")
        fit.init= prc (xvar, dil.x, yvar, dil.y, method=method.init, verbose=verbose) 
        if (!is.null(fit.init$error)) {
            res$error="During initialization, prc fails to run."
            cat(res$error, "\n", file=stderr())
            return(res)
        }
        theta=c(coef(fit.init), sigma.sq=fit.init$sigma.sq)
        theta=c(logc=unname(log(theta["c"])), logd=unname(log(theta["d"])), theta[-match(c("c","d"), names(theta))])
    } else {
        theta=init
        if (length(theta)!=5) stop("init need to include sigma.sq")
    }    
    nu=c(alpha=1,beta=1) # initial value for nu
    if (verbose) {
        cat("\n")
        if (model=="struct") cat("initial nu:", nu, "\n")
        cat("initial theta:", theta, "\n\n")
    }    
    
    # define the grid and compute matrix A for the first time
    u.0=seq(theta["logc"]*0.9, theta["logd"]*1.1, length=2+grid.density*1.2)
    
    # loop through step 2 and 3. 
    # only theta has new.theta because we need to keep both update and previous values to decide convergence
    # we only keep one iteration of support/p/nu
    iterations=0
    histories=list()
    mix.liks.after.step.3=NULL   
    while(TRUE) {
    
        iterations=iterations+1
        if (verbose>=2) cat("=========== Iter "%+%iterations%+%" ==============\n") else cat(formatC(iterations,digits=2,), ") ", sep="")
        
        #########################################################
        # Step 2: update distribution of latent variable either nonparametrically or parametrically
        u=u.0[ u.0<theta["logd"] & u.0>theta["logc"] ] # subset to the grid points in bound
        A=compute.A (theta["logc"], theta["logd"], theta["b"], theta["f"], dil.r, theta["sigma.sq"], u, xvar, yvar)         
        
        #### Nonparametric estimate of distribution
        if (model=="sp") {
            
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
            #mosek_write (sco1, "mosekfiles/yf_ex1_run_1.task", list(scofile="yf_ex1.sco"))
            mix.lik.dual = sum(-log(fit.mosek$sol$itr$xx))
            
            if (verbose>=2) {
                myprint(mix.lik.dual)               
                plot(fit.mosek$sol$itr$xc, sco1$A %*% fit.mosek$sol$itr$xx); abline(0,1) # sanity check
                title(main="Iteration "%+%iterations, outer=T, line=-1)
                empty.plot()
                #print(fit.mosek)
                #print(table(fit.mosek$sol$itr$skc))        
    #            print(table(fit.mosek$sol$itr$skx))
    #            print(summary(fit.mosek$sol$itr$xx))
    #            print(summary(fit.mosek$sol$itr$xc))
    #            print(cbind(fit.mosek$sol$itr$skx, fit.mosek$sol$itr$xx, fit.mosek$sol$itr$slx, fit.mosek$sol$itr$sux))
                if(verbose>=4) print(cbind(fit.mosek$sol$itr$skc, fit.mosek$sol$itr$xc, fit.mosek$sol$itr$slc, fit.mosek$sol$itr$suc)[fit.mosek$sol$itr$skc=="UL",])
            }
            
            #print(sort(fit.mosek$sol$itr$xc, decreasing=TRUE))
        
            # recover primal soln
            support.set.1 = suppressWarnings(which(fit.mosek$sol$itr$skc=="UL")) # UL seems to choose less than ideal set
            if (length(support.set.1)<=1) {
                support.set.1 =order(fit.mosek$sol$itr$suc,decreasing=TRUE) [1:floor(n*.1)] # 10% sample size
                support.set.2 =order(fit.mosek$sol$itr$suc,decreasing=TRUE) [1:floor(n*.2)] # 20% sample size
                support.set.3 =order(fit.mosek$sol$itr$suc,decreasing=TRUE) [1:floor(n*.4)] # 40% sample size
            } else {
                support.set.2 =order(fit.mosek$sol$itr$suc,decreasing=TRUE) [1:min(length(support.set.1)*1.5,n*.8)] # twice as many as support.set.1, but not over 80% sample size
                support.set.3 =order(fit.mosek$sol$itr$suc,decreasing=TRUE) [1:min(length(support.set.1)*2.0,n*.8)] # four times as many as support.set.1, but not over 80% sample size
            }
            # the following line may throw an error when mosek did not converge
            support.set.5 =try(order(fit.mosek$sol$itr$suc,decreasing=TRUE) [1:which(cumsum(sort(fit.mosek$sol$itr$suc,decreasing=TRUE))>=1-1e-6)[1]]) # cumsum of suc just above 1
            if (verbose>=2) cat("support sizes UL|1|2|3: ", suppressWarnings(sum(fit.mosek$sol$itr$skc=="UL")), length(support.set.1), length(support.set.2), length(support.set.3), "\n")
            # Others tried:
            # fit.mosek$sol$itr$xc==n sometimes return empty set
            # a modification of set.1 that includes all with suc greater than any UL suc, does not make a big difference
    #        support.set.2 = which(fit.mosek$sol$itr$suc>0.01) 
    #        support.set.3 = which(fit.mosek$sol$itr$suc>0.001) 
            
            # try three choices of support set and use the one with the best mix.lik.primal
            support.set.list=list()
            if (length(support.set.2)<length(support.set.1)) num.sets=1 else if (length(support.set.3)<length(support.set.2)) num.sets=2 else num.sets=3
            if (!try.additiona.support.sets) num.sets=1
            mix.lik.primals=numeric(num.sets)
            for (support.idx in 1:num.sets) {
                if (verbose>=2) myprint(support.idx)
                support.set = get("support.set."%+%support.idx)            
                while (TRUE) {
                    p.fit=try(lm.fit(x = A[,support.set,drop=FALSE], y = 1/fit.mosek$sol$itr$xx))
                    if (inherits(p.fit, "try-error")) {mix.lik.primal=-Inf; break}
                    p.new=coef(p.fit)
                    if (any(is.na(p.new))) { 
                        #stop("some p.new are NA")
                        p.new[is.na(p.new)]=-1
                    } 
    #               p.new=fit.mosek$sol$itr$suc[support.set] # performance not good
    #                support.1=unique(r.x); p.1=rep(1/5,5) # mix.lik not as good as support
    #                support.1=support[-c(2,4)]; p.1=rep(1/5,5) # mix.lik not as good as support
    #                sum(log(compute.A (theta["logc"],theta["logd"],theta["b"],theta["f"], dil.r, theta["sigma.sq"], support.1, xvar, yvar) %*%p.1))
                    if (verbose>=2) cat("sum(p.new)=", round(sum(p.new),2), ", ", sep="")
                    p.new=p.new/sum(p.new)
                    mix.lik.primal=.Call("compute_mixlik", exp(theta["logc"]),exp(theta["logd"]),theta["b"],theta["f"], dil.r, theta["sigma.sq"], u[support.set], xvar, yvar, p.new)
                    if (verbose>=2) {
                        myprint(mix.lik.primal, digits=6)
                        plot(p.fit$fitted.values, 1/fit.mosek$sol$itr$xx) # making sure we are doing a good job at recovering primal soln
                        plot(0,0, type="n", ylim=range(c(0,p.new)), xlim=range(u[support.set]), xlab="r", ylab="p")#, main="Iteration "%+%iterations)
                        points(u[support.set], p.new, pch=19)
                        for (i in 1:length(p.new)) lines(rep(u[support.set][i],2), c(0,p.new[i]))            
                        if (verbose>=3) {
                            print(cbind(fit.mosek$sol$itr$skc, fit.mosek$sol$itr$xc, fit.mosek$sol$itr$slc, fit.mosek$sol$itr$suc)[order(fit.mosek$sol$itr$suc,decreasing=TRUE)[1:n],])
                            myprint(sum(fit.mosek$sol$itr$suc[support.set]))
                            #myprint(sum(fit.mosek$sol$itr$suc[support.set.5]))
                            print(rbind(u[support.set], p.new))
                        }                            
                    }                    
                    support.set=support.set[p.new>0]
                    if (all(p.new>0)) break 
                } # end while (TRUE)
                mix.lik.primals[support.idx]=mix.lik.primal            
                support.set.list[[support.idx]]=support.set # save this
                if(verbose>=2) {old.par=par(no.readonly = TRUE); par(old.par)} # to start a new page with the existing graphics settings, but the order of mfrow vs mfcol is not saved
            } # end support.idx loop
            mix.lik.primal=max(mix.lik.primals)
            if (verbose>=2) myprint(mix.lik.primal, digits=6)
            if (mix.lik.primal<mix.lik.dual/2) {
                res$error="unable to extract primal soln"
                cat("\n", res$error, "\n", file=stderr())
                break;
            } else {
                # updated support.set has been saved, p needs to be recomputed from that
                support.set = support.set.list[[which.max(mix.lik.primals)]]
                p.fit=lm.fit(x = A[,support.set,drop=FALSE], y = 1/fit.mosek$sol$itr$xx)
                p=coef(p.fit)
                p=p/sum(p)
                support=u[support.set]
            }
                
                
        #### Structural relationship estimate
        } else if (model=="struct"){
            
            support=u
            support.transf = (support-theta["logc"])/(theta["logd"]-theta["logc"])
            # assume gamma distribution
            optim.out = try(optim(
                nu,
                function(nu.fn,...) {
                    p=suppressWarnings(dbeta(support.transf, nu.fn["alpha"], nu.fn["beta"]))
                    p=p/sum(p)
                    -sum(log(A%*%p))
                }, 
                gr=NULL,
                method="BFGS", control = list(trace=0), hessian = F
            ), silent=TRUE)
            if(inherits(optim.out, "try-error")) {
                # try again with Neld-Mead
                optim.out = try(optim(
                    nu,
                    function(nu.fn,...) {
                        p=suppressWarnings(dbeta(support.transf, nu.fn["alpha"], nu.fn["beta"]))
                        p=p/sum(p)
                        -sum(log(A%*%p))
                    }, 
                    gr=NULL,
                    method="Nelder-Mead", control = list(trace=0), hessian = F
                ))   
                if(inherits(optim.out, "try-error")) {
                    res$error="fails to find nu"
                    cat("\n", res$error, "\n", file=stderr())
                    break
                }             
            }
            if (verbose>=4) print(optim.out)
            nu = optim.out$par
            p=dbeta(support.transf, nu["alpha"], nu["beta"])            
            p=p/sum(p)
            if (verbose>=2) {
                cat("opt nu     : ", nu, ";", -optim.out$value, "\n")
            }
            
        }        
        # To summarize, step 2 produces two variables: support and p, and nu if structrual. 
        # For sp model, support is the set of grid points that have non-zero probablity estimate
        # For struct model, support is the grid
    
       
        new.theta=theta        
        max.inner.iter=25
        inner.iter=0        
        debug.1=function(){ # trying to figure out why c_hat from struct has a problem compared to tls estimator
            support.fn=u.0[u.0>new.theta["logc"] & u.0<new.theta["logd"]]
            p.fn=dbeta((support.fn-new.theta["logc"])/(new.theta["logd"]-new.theta["logc"]), nu["alpha"], nu["beta"]); p.fn=p.fn/sum(p.fn)     
            cat(" support: ", length(support.fn), "from", support.fn[1], "to", last(support.fn), ", p[1:2 and last]:", p.fn[1:2], last(p.fn),"\n")
            saved=(log(compute.A (new.theta["logc"], new.theta["logd"], new.theta["b"], new.theta["f"], dil.r, new.theta["sigma.sq"], support.fn, xvar, yvar) %*% p.fn ))
            cat(saved, "\n"); saved        
        }                                
        #saved.old=debug.1()
        while (TRUE) {
            inner.iter=inner.iter+1
            new.theta.0=new.theta # new.theta.0 is used to decide convergence            
    
            #########################################################
            # Step 3a: update curve parameters estimate
            # optim (Nelder-Mead) performs ok, but is too slow
            # optim (BFGS) has trouble with numerical differentiation
            # nlminb performs well, but will return log(c) that is too close to min(support), causing problems for sigmasq if continue to step 3a
            nlminb.out = nlminb(
                new.theta[logdbf],
                function(theta.fn,...) {
                    if (model=="sp") {
                        subset.=support>new.theta["logc"] & support<theta.fn["logd"]                        
                        support.fn=support[subset.]
                        p.fn=p[subset.]
                        p.fn=p.fn/sum(p.fn)
                    } else if (model=="struct") {
                        subset.=u.0>new.theta["logc"] & u.0<theta.fn["logd"]                        
                        support.fn=u.0[subset.]
                        p.fn=dbeta((support.fn-new.theta["logc"])/(theta.fn["logd"]-new.theta["logc"]), nu["alpha"], nu["beta"])            
                        p.fn=p.fn/sum(p.fn)     
                    }
                    # if divide by n, results are slightly different
                    -.Call("compute_mixlik", exp(new.theta["logc"]), exp(theta.fn["logd"]), theta.fn["b"], theta.fn["f"], dil.r, new.theta["sigma.sq"], support.fn, xvar, yvar, p.fn)
                }, 
                gradient=NULL, hessian = F,
                control = list(trace=0)
            )
            new.theta[logdbf] = nlminb.out$par
            if (verbose>=2) cat("opt theta  : ", new.theta, ";", -nlminb.out$objective, "\n")
            if (verbose>=4) print(nlminb.out)            
            #saved.2=debug.1(); plot(saved.2-saved.old); saved.old=saved.2
    
            #########################################################
            # Step 3b: update sigma.sq estimate
            if (model=="sp") {
                subset.=support>new.theta["logc"] & support<new.theta["logd"]                        
                support.3=support[subset.]
                p.3=p[subset.]
                p.3=p.3/sum(p.3)
            } else if (model=="struct") {
                subset.=u.0>new.theta["logc"] & u.0<new.theta["logd"]
                support.3 = u.0[subset.]
                p.3=dbeta((support.3-new.theta["logc"])/(new.theta["logd"]-new.theta["logc"]), nu["alpha"], nu["beta"])            
                p.3=p.3/sum(p.3)
            }
            optimize.out = optimize(
                f = function(sigma.sq.f,...) {
                    -.Call("compute_mixlik", exp(new.theta["logc"]), exp(new.theta["logd"]), new.theta["b"], new.theta["f"], dil.r, sigma.sq.f, support.3, xvar, yvar, p.3)
                },
                interval=c(new.theta["sigma.sq"]/5, new.theta["sigma.sq"]*5), # when the factor is 10, likelihood often drops 
            )
            new.theta["sigma.sq"] = optimize.out$minimum    
            if (verbose>=2) cat("opt sigmasq: ", new.theta, ";", -optimize.out$objective, "\n")
            if (verbose>=4) print(optimize.out)
                    
            
            # stopping rule for step 3
            if (max(abs(1 - new.theta/new.theta.0)) < reltol) {
                break;
            } else if (inner.iter==max.inner.iter) {
                break
            } 
                        
        } # end step 3 loop
        
        mix.lik.after.step.3 = -optimize.out$objective
        mix.liks.after.step.3=c(mix.liks.after.step.3, mix.lik.after.step.3)
        if (model=="sp") histories[[iterations]]=list(support=support, p=p, theta=new.theta, mix.lik.after.step.3=mix.lik.after.step.3)
        if (model=="struct") histories[[iterations]]=list(support=support, p=p, nu=nu, theta=new.theta, mix.lik.after.step.3=mix.lik.after.step.3)
                
        if (verbose>=2) {
            if (model=="sp") cat("support size:", length(support), "\n")            
            cat("\n")            
        } else if (verbose) {
            cat("ll", formatC(mix.lik.after.step.3,format="g",width=7))
            if (model=="sp") cat(", #support", formatC(length(support),format="g",digits=ceiling(log10(n))), sep="")
            if (model=="struct") cat(", nu", nu, sep=" ")
            cat(", theta", new.theta, "\n")                        
        }
        
        # stopping rules
        if (max(abs(1 - new.theta/theta)) < reltol) {if (verbose) cat("converged\n"); break; }        
        theta = new.theta # update theta for the new iteration
        if (iterations>=max.iter) {
            cat("Stop: max iter reached\n"); break;
        } else if (stop.when.dropping & mix.lik.after.step.3<mix.liks.after.step.3[ifelse(length(mix.liks.after.step.3)-1<1,1,length(mix.liks.after.step.3)-1)]) {
            cat("Stop: likelihood starts decreasing.\n"); break;
        } else if (mix.lik.after.step.3 < 0.5*max(mix.liks.after.step.3)) {
            cat("Stop: likelihood decreases to half of the maximum.\n"); break;
        }        
            
    } # end while loop
        
    if (length(mix.liks.after.step.3)==0) {
        return (res)
    } else {
        res$error=NULL # if there are already iterations, we need to remove res$error because this is used to determine overall fit status
    }
    
    # all parameter taken from best.iter in the history
    best.iter=histories[[which.max(mix.liks.after.step.3)]]    
    res$support=best.iter$support
    res$p=best.iter$p
    res$coefficients=best.iter$theta[c("logc",logdbf)]
    res$sigma.sq = best.iter$theta["sigma.sq"]
    res$mixlik = best.iter$mix.lik.after.step.3
    if (model=="struct") {
        res$nu=best.iter$nu
    }
    
    res$A=with(res, compute.A (coefficients["logc"],coefficients["logd"],coefficients["b"],coefficients["f"], dil.r, sigma.sq, support, xvar, yvar) )
    
    res$iterations=iterations    
    if (keep.history) res$histories=histories
    
    #### compute asymptotic variance of theta_hat and xhat and yhat
    
    
        
    # return object
    res
}

# a convenience function
prcstruct=function(xvar, dil.x, yvar, dil.y, 
    grid.density=200, method.init=c("TLS","naive"), reltol=1e-3, max.iter=20, init=NULL,
    keep.history=FALSE, verbose=FALSE) 
{
    method.init=match.arg(method.init)
    prcsp(xvar, dil.x, yvar, dil.y, model="struct",
        grid.density=grid.density, method.init=method.init, reltol=reltol, max.iter=max.iter, init=init,
        keep.history=keep.history, verbose=verbose)
}

# this should match the mixture likehood obtained directly from mosek
mixlik <- function(object, ...) UseMethod("mixlik") 
mixlik.prc=function(object, ...) {    
    #out.1=with(object, sum(log(compute.A (coefficients["logc"],coefficients["logd"],coefficients["b"],coefficients["f"], dilution.ratio, sigma.sq, support, xvar, yvar) %*%p)))
    out.2=with(object, .Call("compute_mixlik", exp(coefficients["logc"]),exp(coefficients["logd"]),coefficients["b"],coefficients["f"], dilution.ratio, sigma.sq, support, xvar, yvar, p))
    out.2
}

# populate A
compute.A = function(logc,logd,b,f, dil.r, sigma.sq, support, xvar, yvar) {
    c=exp(logc); d=exp(logd)
    .Call("compute_A", c,d,b,f, dil.r, sigma.sq, support, xvar, yvar)
    # .Call is a lot faster than the R implementation below
#    m.f.1=function(c,d,b,f,r,x,y,dil.r)  (y - four_pl_prc(c,d,b,f, r, dil.r))^2  +  (x - r)^2 # this is about four times as fast as m.f
#    A=matrix(NA,n,K)
#    for (i in 1:n) {
#        for (k in 1:K) {
#            m_ik = m.f.1(c, d, b, f, u[k], xvar[i], yvar[i], dil.r) 
#            A[i,k] = 1/fit$sigma.sq * exp(-m_ik/fit$sigma.sq/2)
#        }
#    }
}


compute.m = function(logc,logd,b,f, dil.r, sigma.sq, support, xvar, yvar) {
    c=exp(logc); d=exp(logd)
    n=length(xvar)
    m.f.1=function(c,d,b,f,r,x,y,dil.r)  (y - four_pl_prc(c,d,b,f, r, dil.r))^2  +  (x - r)^2 # this is about four times as fast as m.f
    K=length(support)
    m=matrix(NA,n,K)
    for (i in 1:n) {
        for (k in 1:K) {
            m[i,k] = m.f.1(c, d, b, f, support[k], xvar[i], yvar[i], dil.r) 
        }
    }
    m
}


###### some code used to be at the begining of step 3
#        # initialize new.theta with an approximation to the mixture likelihood, this seems quaint a
#        A=compute.A (theta["logc"], theta["logd"], theta["b"], theta["f"], dil.r, theta["sigma.sq"], support, xvar, yvar) 
#        best.support=support[apply(A, 1, which.max)]
#        s.sqrt.k.rvar = four_pl_prc (exp(theta["logc"]), exp(theta["logd"]), theta["b"], theta["f"], best.support, dil.r^.5) # an easy to make mistake is to use k*.5 instead of k^.5
#        dat.stacked=data.frame(readout=c(unname(xvar),unname(yvar)), x=rep(s.sqrt.k.rvar,2), k=rep(c(dil.r^(-1/2),dil.r^(1/2)),each=n))
#        formula.gnls = as.formula(  "(readout) ~ log(exp(logc)+(exp(logd)-exp(logc))/(1+k^b*(((exp(logd)-exp(logc))/(exp(x)-exp(logc)))^(1/f)-1))^f)"  ) 
#        fit.1=try(gnls(formula.gnls, data=dat.stacked, start=theta[c("c",logdbf)], # key to take only cdbf from new.theta as start for gnls
#            control=gnlsControl(
#            nlsTol=1e-1,  # nlsTol seems to be important, if set to 0.01, then often does not converge
#            tolerance=1e-4, 
#            # msTol=1e-1, minScale=1e-1, .relStep=1e-7,
#            returnObject=TRUE, # allow the return of the fit when the max iter is reached
#            maxIter=5000, nlsMaxIter=50, opt="nlminb", msVerbose=T)), silent=FALSE
#        )       
#        new.theta=coef(fit.1)
#        new.theta=c(new.theta, theta["sigma.sq"]) # has to be on two lines, otherwise names are lost
#        if (any(support<new.theta["logc"]) | any(support>new.theta["logd"])) {
#            #some support is outside new.theta, revert new.theta to theta
#            new.theta=theta
#            if (verbose>=2) cat("initial theta for step 3:", new.theta, "\n")        
#        } else {
#            if (verbose>=2) cat("initial theta for step 3:", new.theta, "\n")                    
#        }
        
