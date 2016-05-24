krm.most = function(formula, data, regression.type=c("logistic","linear"), 
    kern.type=c("rbf","mi","mm","prop"), n.rho=10, range.rho=0.99, 
    n.mc=2000, 
    seq.file.name=NULL, formula.kern=NULL, seq.start=NULL, seq.end=NULL,
    inference.method=c("parametric.bootstrap", "perturbation", "Davies"),
    verbose=FALSE) 
{
    begin = Sys.time()
    
    # match input parameters
    if (length(regression.type)!=1) stop("regression.type has to be specified")
    regression.type <- match.arg(regression.type)    
    kern.type <- match.arg(kern.type)    
    inference.method <- match.arg(inference.method)
    if (inference.method=="Davies") {
        if (kern.type=="mi") {
            cat("Profile HMM MI kernel is currently not implemented with Davies\n")
            return (NA) # the way rho is defined for Davies does not apply to MI kernel
        }
        n.rho=500 # according to Davies paper
    }
    
    # check input parameters
    if (!is.null(seq.file.name) & !kern.type %in% c("mi","mm","prop")) stop("choose a kernel from mi, mm, prop for sequence kernels")
    if (is.null(formula.kern) & is.null(seq.file.name)) stop("either formula.kern or seq.file.name has to be supplied")
    if (!is.null(seq.file.name)) {stopifnot (is.character(seq.file.name)); if(verbose) myprint(seq.file.name) }
    if (!is.null(formula.kern)) {stopifnot (is(formula.kern,"formula")); if(verbose) {cat("formula.kern: "); print(formula.kern)}}
    if (!is.null(seq.start) & !is.null(seq.end) & !is.null(seq.file.name)) {stop("seq.start and seq.end are only supported when sequences are specified through form.kern")}
    
    # save rng state before set.seed in order to restore before exiting this function
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {        
        set.seed(1)
        save.seed <- get(".Random.seed", .GlobalEnv)
    }                        
    
    # compute kernel for rho=1
    if (kern.type %in% c("mi","mm","prop")) {
        if (!is.null(seq.file.name)) {
            K.1 = getSeqKernel(seq.file.name, kern.type, tau = 1, call.C = T)
        } else {
            tmp=as.character(formula.kern)
            if (length(tmp)!=2) stop("formula.kern needs to be in the form ~ seq.col")
            seq.list = as.list(data[[tmp[2]]])
            K.1 = getSeqKernel(seq.list, kern.type, tau = 1, call.C = T, seq.start=seq.start, seq.end=seq.end)
        }
        
    } else {
        Z=model.matrix(formula.kern, data)
        if (colnames(Z)[1]=="(Intercept)") Z=Z[,-1,drop=FALSE]
        K.1 = kyotil::getK(Z, kernel=kern.type, para=1) 
    }    
    
    # define rhos
    l2dist = -log(K.1)
    diag(l2dist)=NA
    if (inference.method!="Davies") {
        median.log.K = median(l2dist, na.rm=T)
        if (median.log.K==0) {print("l2dist has mean and median 0"); return (rep(NA,4))}    
        rho.alpha=(1-range.rho)/2
        rhos=sqrt(-median.log.K/log(seq(rho.alpha,1-rho.alpha, length=n.rho))) # try to space rho uniformly on the correlation scale
        
    } else {
        # implementing Liu et al 2008
        maxk = max(l2dist, na.rm=T)
        mink = min(l2dist, na.rm=T)
        if (mink==0) {
            l2dist[l2dist==0]=NA
            mink = min(l2dist, na.rm=T)
        }
        rhos=sqrt(seq(mink/10, maxk*10, length=n.rho))
    }
    if (verbose) {
        myprint(summary(c(l2dist)))
        myprint(summary(c(rhos)))
        t.1=Sys.time()
        cat("rhos defined. time passed: "); print(t.1-begin)
    }
    
    # terms used in parametric boostrap and perturbation
    X=model.matrix(formula, data)
    X.fr=model.frame(formula, data)[-1]
    y=model.frame(formula, data)[[1]]
    if (verbose==2) {cat("str(y):\n"); str(y)}
    n=nrow(data)
    if (regression.type=="logistic") {    
        fit=glm(formula, data, family="binomial") # null model fit
        # terms derived from the model and don't depend on kernel, many are not needed for parametric bootstrap, but computed anyway
        beta.h=coef(fit)
        mu.h=drop(expit(X %*% beta.h))
        D.h = c(mu.h*(1-mu.h)) 
        XD.5 <- X * D.h^.5
        . <- eigen(crossprod(XD.5))
        V.beta.h <- crossprod(t(.$vectors) * .$values^(-0.5))  # solve(t(X) %*% D.h %*% X)
        V.eta.h  <- crossprod(t(X %*% .$vectors) * .$values^(-0.5)) # X %*% V.beta.h %*% t(X)
        V.mu.h   <- DXD(D.h,V.eta.h,D.h) # D.h %*% V.eta.h %*% (D.h)
        P.h1=diag(D.h) - V.mu.h # this is the variance of Y - mu^hat to be used in simulation mvn
        extra.kurtosis = (mu.h*(1-mu.h)^4 + (1-mu.h)*mu.h^4) - 3 * D.h**2 # needed for mean estimate
        A.h=diag(n) - DXD(D.h,V.eta.h,rep(1,n))# D.h %*% V.eta.h, needed for variance estimate        
            
    } else if (regression.type=="linear") {    
        fit=glm(formula, data, family="gaussian") # null model fit
        # terms derived from the model and don't depend on kernel, many are not needed for parametric bootstrap, but computed anyway
        beta.h=coef(fit)
        mu.h=drop(X %*% beta.h)
        noise.sd = summary(fit)$dispersion ** .5
        P = diag(n) - X %*% solve(crossprod(X)) %*% t(X)        
    
    }     
    resid.=y-mu.h
    if(verbose) {myprint(dim(X)); myprint(beta.h)}
    
    
    # do inference
    Q.rho.stats = sapply(rhos, simplify = "array", function(rho) {
        if (verbose==2) myprint(rho)
        
        # compute K_rho
        K = K.1 ^ (rho^-2)
    
        # parametric bootstrap-based inference
        if (inference.method=="parametric.bootstrap"){        
            test.stats.obs = krm.score.test (formula, data, K, regression.type) 
            
            test.stats.rep = sapply(1:n.mc, simplify = "array", function(j){
                if (verbose>=2) myprint(j)
                set.seed(j+1e6)                         
                if (regression.type=="logistic") {
                    new.y <- rbern(n, mu.h)
                } else if (regression.type=="linear") {
                    new.y <- rnorm(n, mean=mu.h, sd=noise.sd) 
                }
                new.dat = cbind(new.y, X.fr)
                names(new.dat)[1] = as.character(formula)[2]
                krm.score.test (formula, new.dat, K, regression.type, verbose=verbose) 
            })
            cbind(test.stats.obs,test.stats.rep)        
        
        # perturbation-based inference 
        } else if (inference.method=="perturbation") {
            requireNamespace("MASS")
            if (regression.type=="logistic") {            
                m1=tr(P.h1 %*% K)
                str(K)
                str(A.h)
                W.h=crossprod(A.h,symprod(K,A.h)) # t(A.h) %*% K %*% A.h
                v2 = varQ (W.h, variance=D.h, extra.kurtosis=extra.kurtosis, do.C=TRUE)
                test.stats.obs = perturbation.test(resid., K, m1, v2)  
                if (verbose) {
                    test.stats.obs.2 = krm.score.test (formula, data, K, regression.type, verbose=FALSE)           
                    stopifnot(all(test.stats.obs.2==test.stats.obs, na.rm=TRUE)) # sanity check
                }
    
                test.stats.rep = sapply(1:n.mc, simplify = "array", function(j){
                    set.seed(j+1e6) 
                    new.resid <- MASS::mvrnorm(1,mu=rep(0,n),Sigma=P.h1) 
                    perturbation.test(new.resid, K, m1, v2)                        
                })
            
            } else if (regression.type=="linear") {
                PK <- P %*% K            
                m=tr(PK)        
                V.Q.norm = 2*tr(crossprod(PK)) # 2*tr( P %*% K %*% P %*% K )
                test.stats.obs = perturbation.test(resid./noise.sd, K, m, V.Q.norm)  
                if (verbose) {
                    test.stats.obs.2 = krm.score.test (formula, data, K, regression.type)
                    stopifnot(all(test.stats.obs.2==test.stats.obs, na.rm=TRUE)) # sanity check
                }
                               
                test.stats.rep = sapply(1:n.mc, simplify = "array", function(j){
                    set.seed(j+1e6) 
                    new.resid <- MASS::mvrnorm(1,mu=rep(0,n),Sigma=P) 
                    perturbation.test(new.resid, K, m, V.Q.norm)                        
                })
                
            }
            cbind(test.stats.obs,test.stats.rep)        
            
        } else if (inference.method=="Davies") {
            if (regression.type=="logistic") {            
                Ph1K <- P.h1 %*% K
                m1=tr(Ph1K)
                V.Q.norm.h = 2*tr(crossprod(Ph1K)) # 2*tr( P.h1 %*% K %*% P.h1 %*% K )
                Q = txSy(resid.,K,resid.) # drop(t(new.resid) %*% K %*% new.resid)    
                (Q-m1)/sqrt(V.Q.norm.h)
                
            } else stop("Davies is logistic regression only")
        
        } else stop ("something wrong in krm.most")
                
        
    }) 
    ## if inference.method is not Davies, Q.rho.stats is 3-dimensional array: # of test stats   by   1+n.mc   by   # of rhos
    ## if inference.method is Davies, Q.rho.stats is vector of length n.rho
    assign(".Random.seed", save.seed, .GlobalEnv) # restore rng state 
    if (verbose==2) {cat("str(Q.rho.stats):\n"); str(Q.rho.stats)}
        
    end=Sys.time()
    cat("Total time passed: "); print(end-begin)
    
    if (inference.method %in% c("parametric.bootstrap","perturbation")) {
        # maximize over # of rhos
        p.chi.sup.1 <- drop(apply(Q.rho.stats["chiI",,,drop=FALSE],1:2,max))
        p.chi.sup.2 <- drop(apply(Q.rho.stats["chiII",,,drop=FALSE],1:2,max))
        Q.norm.sup.1 <- drop(apply(Q.rho.stats["normI",,,drop=FALSE],1:2,max))
        Q.norm.sup.2 <- drop(apply(Q.rho.stats["normII",,,drop=FALSE],1:2,max))    
        if (verbose==2) {cat("str(p.chi.sup.1):\n"); str(p.chi.sup.1)}    
            
        p.values = c(
            "chiI"=mean(p.chi.sup.1[-1] > p.chi.sup.1[1]), 
            "chiII"=mean(p.chi.sup.2[-1] > p.chi.sup.2[1]), 
            "normI"=mean(Q.norm.sup.1[-1] > Q.norm.sup.1[1]), 
            "normII"=mean(Q.norm.sup.2[-1] > Q.norm.sup.2[1])
        )
    } else if (inference.method=="Davies") {
        M = max(Q.rho.stats)
        W = sum(abs(diff(Q.rho.stats)))
        p.values = pnorm(-M) + W*exp(-0.5*(M^2))/sqrt(8*pi)
    }
    
    out=list()
    class(out)=c("krm", class(out))
    out$p.values=p.values
    #out$Q.rho.stats=Q.rho.stats
    out
}


#
perturbation.test = function (new.resid, K, m, v) {    
    Q = txSy(new.resid,K,new.resid) # drop(t(new.resid) %*% K %*% new.resid)    
    s=v/(2*m); k=2*m^2/v
    c("chiI"=pchisq(Q/s, df=k), "chiII"=NA, "normI"=pnorm((Q-m)/sqrt(v)), "normII"=NA)    
}
