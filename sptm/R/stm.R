# if a clinical covariate has missingness, program fails, to be fixed
# when verbose is >= 3, it is debug mode
stm = function (formula, dat, strata.formula, phase2.ind=NULL, imputation.formula=NULL, 
    family=c("PH","PO","P2"), ee=c("fine2","fine1","kong"), var.est.type=c("1","2"),
    t0, init=NULL, maxit=1000,
    intermediate=FALSE, verbose=FALSE, show.time.elapsed=TRUE)
{
    
    begin = Sys.time()
    
    ee <- match.arg(ee)  
    if (length(family)>1) stop("The family argument has to be specified, choose from PH, PO and P2")  
    family <- match.arg(family)    
    var.est.type <- match.arg(var.est.type)    
        
    #### infer phase2.ind
    if (is.null(phase2.ind)) {
        vars.in.formula=intersect(all.vars(formula), attr(terms(formula),"term.labels")) 
        phase2.ind=complete.cases(dat[,vars.in.formula])
    } else {
        stopifnot(nrow(dat)==length(phase2.ind))
        stopifnot(is.logical(phase2.ind))
    }
    
    # fit cox model to get init values
    if (is.null(init)) {
        cox.fit = svycoxph(formula, design=twophase(id=list(~1,~1),strata=list(NULL,strata.formula),subset=phase2.ind,data=dat))
        init=c(log(0.032*t0),coef(cox.fit))
    }        
    
    if (verbose>=3) print("before calling stm.internal")
    if (verbose>=3) myprint(init)    
    # this provides a init value for cailbration fit, especially the first parameter. this improves convergence when doing calibration
    res = stm.internal(formula, dat, strata.formula, phase2.ind=phase2.ind, 
        family=family, ee=ee, var.est.type=var.est.type,
        t0=t0, init=init, maxit=maxit,
        intermediate=intermediate, verbose=verbose, show.time.elapsed=FALSE)
    if (verbose>=3) print("after calling stm.internal")
    
    if (is.null(imputation.formula)) {
        # no more to be done
    } else {
        
        init=res[,"est"]
        
        ####################################################################################################
        ### calibration is done in three steps
        
        #########################################
        #### step 1 (a)
        if (verbose) cat("------------------------------------------------------------------------------------------")
        if (verbose) cat("\nStep 1(a): Predict missing phase 2 covariates using imputation model\n")
#        # use lm to impute
#        fit.step1=lm(imputation.formula, dat)
        # use svyglm to impute
        dstrat<-twophase(id=list(~1,~1),strata=list(NULL,strata.formula),subset=phase2.ind,data=dat) 
        fit.step1 = svyglm(imputation.formula, design=dstrat)
        # make prediction for all observations, not just the ones with missing covariates
        predicted = predict(fit.step1,type="response",newdata=dat,se=F)
        # the left hand side in the imputation model might be a variable name, e.g. s, or a transformation, e.g. logit(s)
        lhs=as.character(imputation.formula)[2]
        dat.step1 = dat 
        if (contain(lhs, "(")) {
            tmp=strsplit(lhs,"[\\()]")[[1]]
            transf = tmp[1]
            lhs = tmp[2]
            if (transf=="logit") transf.f=expit else stop ("transformation not supported")
            dat.step1[,lhs] <- transf.f(predicted)
        } else {
            dat.step1[,lhs] <- predicted
        }
        t.step1a=Sys.time()
        if(verbose) {   
            cat("Time spent in step 1(a): "%+%format(t.step1a-begin)%+%"\n")
        }
        
        
        #########################################
        #### step 1 (b) 
        if (verbose) cat("\n------------------------------------------------------------------------------------------")
        if (verbose) cat("\nStep 1(b): Fit augmented dataset with interest model to get efficient influence function\n")
        # using cox estimate as init improves convergence speed just a little
        # variance is not computed when return.eif is TRUE
        init1 = c(init[1], coef(coxph(formula, dat.step1))) # this has better convergence property than init1=init
        p = length(init1) # number of regression coefficients, including intercept, defined here because p is used outside while loop and there is a break at the end of this block
        fit.step2 = stm.internal(formula, dat.step1, strata.formula, phase2.ind=rep(TRUE,nrow(dat.step1)), family=family, 
            t0=t0, init=init1, maxit=maxit,
            return.eif=TRUE, ee=ee, intermediate=intermediate, show.time.elapsed=FALSE, verbose=verbose)
        r=attr(fit.step2, "eif")
        if(is.null(r)) {
            cat ("\nStep 1(b), fitting augmented dataset, fails to converge.\n")
            res=NULL; 
        } else {
            
            t.step1b=Sys.time()
            if(verbose) cat("Time spent in step 1(b): "%+%format(t.step1b-t.step1a)%+%"\n")
            
            
            #########################################
            #### step 2
            if (verbose) cat("\n------------------------------------------------------------------------------------------")
            if (verbose) cat("\nStep 2: Doing weight calibration\n")
            
            N=nrow(dat)
            n=sum(phase2.ind)
            II = rep (1:N, N)
            JJ = rep (1:N, each=N)        
            ii = rep (1:n, n)      # 1 2 1 2
            jj = rep (1:n, each=n) # 1 1 2 2
            
            # compute weights from strata.formula and phase2.ind
            sampling.p = ave(phase2.ind, model.frame(strata.formula,data=dat), FUN=sum)/
                         ave(phase2.ind, model.frame(strata.formula,data=dat), FUN=length)
            sampling.p.ph2 = sampling.p[phase2.ind] # sampling.p is vector of length N, sampling.p.ph2 is vector of length n
            weights.ph2=1/sampling.p.ph2[ii]/sampling.p.ph2[jj] # dim n^2 by 1
            
            # solve the raking estimating equations to estimate lambda
            xi=phase2.ind[II] & phase2.ind[JJ]        
            r.ph2 = r[xi,] # dim n^2 by p
            r.colsums = colSums(r)
            raking=function(lambda) {
                est.func = drop(weights.ph2 * exp(- r.ph2 %*% lambda)) * r.ph2
                sum( (colSums(est.func) - r.colsums) ** 2 )
            }
            abstol=0.01
            raking.optim=optim(rep(0,p), raking, control=list(abstol=abstol, reltol=1e-100, maxit=10000))        
            if (raking.optim$value>=abstol) {
                cat ("Raking fails to converge. \n")
                str(raking.optim)
                res=NULL
            } else {
                
                lambda = raking.optim$par  
                if (verbose) cat ("calibration coefficient - lambda: ", lambda, "\n")
                       
        #        # alternatively, one may find lambda by approximation formula, faster but sometimes not as accurate
        #        lambda=c(solve(t(r) %*% r) %*% colSums((rho-1)*r))
                
        #        # examine lambda
        #        r.bar = colMeans(r)
        #        b1 = solve((t(r) %*% r)/N/N)
        #        lambda.hat.expansion = b1 %*% colMeans(xi/selection.prob*r - rep(r.bar, each=N*N))
        #        cbind(lambda, lambda.hat.expansion)
        #        summary(1-r%*%lambda.hat.expansion)
        #        summary(exp(- r %*% lambda)-1 - -(r%*%lambda))
        #        summary(exp(- r %*% lambda.hat.expansion)-1 - -(r%*%lambda.hat.expansion))
            
                t.step2=Sys.time()
                if(verbose) cat("Time spent in step 2: "%+%format(t.step2-t.step1b)%+%"\n")
                        
            
                #########################################
                #### step 3
                if (verbose) cat("\n------------------------------------------------------------------------------------------")
                if (verbose) cat("\nStep 3: Estimate regression coefficients\n")
                res = stm.internal(formula, dat, strata.formula, phase2.ind=phase2.ind, family=family,
                    t0=t0, init=init, maxit=maxit,
                    return.eif=FALSE, ee=ee, intermediate=intermediate, var.est.type=var.est.type,
                    verbose=verbose, show.time.elapsed=FALSE,
                    calibration=list(lambda=lambda, r=r))
            
                t.step3=Sys.time()
                if(verbose) cat("Time spent in step 3: "%+%format(t.step3-t.step2)%+%"\n")
                
            } # end if is.null(r)
    
        } # end if raking.optim$value>=abstol
                        
    } # end     if (is.null(imputation.formula)) 
    
    end=Sys.time()
    if(show.time.elapsed) cat("\nTotal time passed: "%+%format(end-begin)%+%"\n")
    
    if(is.null(res)) res = matrix(NA, nrow=p, ncol=2)
    
    return (res)
}


# the first element of init should be log(baseline.hazard)
stm.internal = function (formula, dat, strata.formula, phase2.ind,
    family, ee, var.est.type,
    t0, init=NULL, maxit=1000,
    intermediate=FALSE, return.eif=FALSE, verbose=FALSE, show.time.elapsed=TRUE,
    calibration=NULL)
{
    
    begin = Sys.time()
    
    # infer Xcol and dcol. Xcol is the name of the column for the shorter of the failure and censoring time, dcol is the column name of the disease outcome
    tmp.s = strsplit(as.character(formula)[2], " *\\( *| *, *| *) *")[[1]]
    if (tmp.s[1]!="Surv") stop ("something is wrong with the formula for the interest model")
    Xcol=tmp.s[2]
    dcol=tmp.s[3]
    
    dat.phase2=dat[phase2.ind,]     
    N=nrow(dat) # N is the number of subjects in the full cohort    
    N.1 = sum(dat[[dcol]])    
    N.0 = N - N.1
    n=sum(phase2.ind) # n is the number of subjects in the subcohort
    n.1 = sum(dat[[dcol]] & phase2.ind)    
    ii = rep (1:n, n)      # 1 2 1 2
    jj = rep (1:n, each=n) # 1 1 2 2
    
    #### compute weights from strata.formula and phase2.ind
    sampling.p = ave(phase2.ind, model.frame(strata.formula,data=dat), FUN=sum)/
                 ave(phase2.ind, model.frame(strata.formula,data=dat), FUN=length)
    sampling.p.ph2 = sampling.p[phase2.ind] # sampling.p is vector of length N, sampling.p.ph2 is vector of length n
    weights.ph2=1/sampling.p.ph2[ii]/sampling.p.ph2[jj]
    
    #### estimate survival function of censoring variable
    # estimate censoring function with the full cohort, because (X,d) are available for all
    censoring.s = survfit(as.formula("Surv("%+%Xcol%+%",1-"%+%dcol%+%")~1"), data=dat)    
    idx=sapply(1:n, function(j) sum(censoring.s$time-dat.phase2[[Xcol]][j] <= 1e-10) ) # if we use censoring.s$time <= dat.phase2[[Xcol]][j], we wil run into numerical problem
    # choose one of three
    # following Cheng et al. and Fine et al. and use survival function estimates, which is Pr(C>x)
    G.hat = (censoring.s$surv)[idx] 
#    # to get Pr(C>=x) by moving idx down by 1
#    G.hat = ifelse(idx>1, censoring.s$surv[idx-1], 1) 
#    # true survival function 
#    #G.hat = ifelse(dat.phase2$X<3, 1-dat.phase2$X/15, 0.8) # for sim.fong
#    G.hat = ifelse(dat.phase2$X<.185, 1-dat.phase2$X/.185, 0) # for sim.kong
    
    # G.hat=Pr(C>x) and is 0 at the end of study if there is no random censoring
    A.ij = ifelse(dat.phase2[[Xcol]][jj]<=t0, dat.phase2[[dcol]][jj]*(dat.phase2[[Xcol]][ii]>=dat.phase2[[Xcol]][jj])/G.hat[jj]**2, 0)
    
    #### create building blocks for estimating functions
    # remove intercept column
    Z=model.matrix(formula, dat.phase2)[,-1,drop=F] 
    p=ncol(Z)+1 # Z does not have intercept
    if (n!=nrow(Z)) stop ("n!=nrow(Z) " %+% n %+% " != " %+% nrow(Z))            
    Z.ij = Z[ii,,drop=F]-Z[jj,,drop=F] # n**2 by p-1
    Z.ii = Z[ii,,drop=F] # n**2 by p-1
    Z.jj = Z[jj,,drop=F] # n**2 by p-1
    
    
    #### When stm.cal() calls stm.internal(), it will pass a list as calibration
    raking=!is.null(calibration)
    if (raking) {
        stopifnot(names(calibration)==c("lambda","r"))
        r=calibration$r
        II = rep (1:N, N)
        JJ = rep (1:N, each=N)
        xi=phase2.ind[II] & phase2.ind[JJ]
        r.ph2 = r[xi,] # dim n^2 by p
        weights.ph2.unmod = weights.ph2
        weights.ph2 = weights.ph2 * drop(exp(- r.ph2 %*% calibration$lambda))
        rho=xi
        rho[xi]=weights.ph2.unmod    
        b2.0 = ((rho-1)*r) %*% solve(crossprod(r)) # used in computing h.rho.minus.1, should not compute b2.0 %*% t(r.ph2) since that matrix has a dimension of N^2 by n^2 and cannot fit in memory
    }
    
    #########################################################
    # print debug info
    if (verbose>=2) {
        myprint(N, N.1, n, n.1, n*n, p)
        print(table(sampling.p.ph2))
        if (!raking) print(table(weights.ph2))
        myprint(summary(A.ij))
        if (raking) cat("Raking sampling weights --------\n")
    }
    if (verbose>=4) {
        myprint(censoring.s$time)
        myprint(censoring.s$surv)
        print(cbind(dat.phase2[[Xcol]][jj], G.hat[jj], A.ij, jj))
    }

    # compute any terms that do not depend on theta outside calc to improve performance
    # try to avoid memory operations associated with cbind in calc() by creating B and eta.dot too
    B = cbind(1, Z.ij) 
    ii.neq.jj = ii != jj
    ii.eq.jj = ii == jj
    #eta.dot = matrix(0,n^2,p, dimnames=list(NULL, c("(Intercept)",colnames(Z))))

    # compute estimate, efficient influence function, or variance, when type = 1,2,3, respectively
    # this is implemented as an inner function because it saves us from computing Z.ij etc repeatedly during optim
    # type is default to 1 so that optim only needs to supply theta
    calc = function (theta, type=1) {
        
#        if (verbose>=3) print(theta)
#        if (verbose>=3) str(Z)
        alpha=theta[1]; beta=matrix(theta[-1], ncol=1)
        alpha.exp=exp(alpha)                     
        Z.beta = drop(exp(Z%*%beta))
        Z.beta.ii = Z.beta[ii]
        Z.beta.jj = Z.beta[jj]
        Z.beta.ij = Z.beta.ii/Z.beta.jj
        
        # colnames of eta.dot get passed to the final results
        if (family=="PH") {
            a1 = 1/(1 + Z.beta.ij)
            a2 = alpha.exp*Z.beta.jj/ a1
            a3 = Z.ij / c(1+exp(-Z.ij%*%beta))
            eta.ij = a1*(1-exp(-a2))  
            eta.dot = a1* cbind( "(Intercept)"=a2*exp(-a2),  exp(-a2)*(Z.ii*c(exp(alpha+Z.ii%*%beta)) + Z.jj*c(exp(alpha+Z.jj%*%beta)) + a3) - a3 ) 
        
        } else if (family=="PO") {
            ai = 1 + alpha.exp*Z.beta.ii
            aj = 1 + alpha.exp*Z.beta.jj
            a3 = 1 - Z.beta.ij
            eta.ij = 1/a3*(1-1/aj) - (1-a3)/a3/a3*(log(aj)-log(ai))
            eta.dot = cbind( "(Intercept)"=(aj-1)/ai/aj/aj,  
                        Z.jj*((aj-1)/aj/aj/a3) + (ai-1)/a3/a3*(Z.ii*((1-a3)/ai) + (Z.ii-2*Z.jj)/aj) - Z.ij*((1-a3)*(2-a3)/a3/a3/a3*(log(aj)-log(ai))) )
            eta.ij[ii.eq.jj]=0
            eta.dot[ii.eq.jj,]=0
            # take care of the special case
            aj.1 = aj[drop(Z.ij%*%beta==0)]
            eta.ij[drop(Z.ij%*%beta==0)]  =-0.5*(1/aj.1^2-1)
            eta.dot[drop(Z.ij%*%beta==0),]=cbind("(Intercept)"=1/aj.1^3*(aj.1-1), Z.jj[drop(Z.ij%*%beta==0),]/aj.1^3*(aj.1-1))
                        
        } else if (family=="P2") {
            ai = 1 + 2*c( exp(alpha+Z.ii%*%beta) )
            aj = 1 + 2*c( exp(alpha+Z.jj%*%beta) )
            a3 = 1 - c( exp(Z.ij%*%beta) )        
            eta.ij = 1/a3-1/a3*sqrt(ai/aj)
            eta.dot = cbind( "(Intercept)"=(aj-1)/2/ai^.5/aj^1.5,  
                      Z.ij*(1-a3)/a3^2 - Z.ij*(1-a3)/a3^2*ai^.5/aj^.5 - Z.ii/2/a3*(ai-1)/ai^.5/aj^.5 + Z.jj*(aj-1)/2*ai^.5/a3/aj^1.5 )
            eta.ij[ii.eq.jj]=0
            eta.dot[ii.eq.jj,]=0            
            # take care of the special case
            aj.1 = aj[drop(Z.ij%*%beta==0)]
            eta.ij[drop(Z.ij%*%beta==0)]  =(1-1/aj.1)
            eta.dot[drop(Z.ij%*%beta==0),]=cbind("(Intercept)"=1/aj.1^2*(aj.1-1), Z.jj[drop(Z.ij%*%beta==0),]/aj.1^2*(aj.1-1))
        } 
        
        if (ee=="fine1") {        
            B[,1] = eta.dot[,1] 
        } else if (ee=="kong") {        
            B[,] = eta.dot / eta.ij / (1-eta.ij)
        } # no need to update anything for fine2
        
        E = ii.neq.jj * (A.ij-eta.ij) * B # E is U minus the weights.ph2, used for computing V.ph2
        U = weights.ph2 * E # U is a matrix of n^2 x p
        gr = - crossprod(weights.ph2 * ii.neq.jj * B, eta.dot)
                
        if (type==1) { 
            #### compute estimate
            out = sum (colSums(U)**2)
        
        } else if (type==2) {                                                            
            #### compute analytical variance
            
            #h.rho.minus.1=NULL # only needed for some intermediate checking, not important, comment or remove after using
            
            if (!raking) {     
                
                if (var.est.type==1) {
                    # the first variance estimator, directly estimate it in one step
                    V1=matrix(0,p,p)
                    tmp = .C("V1", as.integer(n), as.integer(p), .as.double((U)), .as.double(V1))
                    V1 = matrix(tmp[[4]], nrow=p)
                    
                } else if (var.est.type==2) {
                    # the second variance estimator, based on Kong et al 2004, break into phase 1 and phase 2
                    V.ph1=matrix(0,p,p)
                    tmp = .C("V1_ph1", as.integer(n), as.integer(p), .as.double((E)), .as.double(V.ph1), .as.double(sampling.p.ph2))
                    V.ph1 = matrix(tmp[[4]], nrow=p) # /N^3 * n^3
                    
                    phi=matrix(0,n,p)
                    tmp = .C("get_phi", as.integer(n), as.integer(p), .as.double((U)), .as.double(phi))
                    phi = matrix(tmp[[4]], nrow=n) * sampling.p.ph2
                    V.ph2 = cov(phi[dat.phase2[[dcol]]==0,]) * (N/n - 1) * N.0 # * (1/n - 1/N) * N.0 / N^2 * N^3
                    
                    V1 = V.ph1 + V.ph2
                }    
                
                V2=matrix(0,p,p)
                tmp = .C("V2", as.integer(n), as.integer(p), .as.double(t(weights.ph2*B*(ii!=jj)*A.ij)), .as.double(V2), 
                    .as.double(dat.phase2[[Xcol]][jj]), as.integer(N), as.integer(dat[[dcol]]), .as.double(dat[[Xcol]]))
                V2 = -4 * matrix(tmp[[4]], nrow=p)
                
                if (verbose==3) {
                    myprint(diag(V1))
                    myprint(diag(V2))
                }
                    
                
            } else { # raking
                        
                e.rho.ph2 = weights.ph2.unmod * E
                h.rho.minus.1 = b2.0 %*% crossprod(r.ph2, e.rho.ph2)
                
                # pad up e.rho.ph2
                e.rho=matrix(0,nrow=N*N,ncol=ncol(e.rho.ph2))
                e.rho[xi,]=e.rho.ph2
                                        
                V1=matrix(0,p,p)
                tmp = .C("V1", as.integer(N), as.integer(p), .as.double((e.rho-h.rho.minus.1)), .as.double(V1)) # the right one
                #tmp = .C("V1", as.integer(n), as.integer(p), .as.double((U)), .as.double(V1)) # U
                #tmp = .C("V1", as.integer(N), as.integer(p), .as.double((h.rho.minus.1)), .as.double(V1)) # check h.rho.minus.1
                V1 = matrix(tmp[[4]], nrow=p)
    
                V2=matrix(0,p,p)
                tmp = .C("V2", as.integer(n), as.integer(p), .as.double(t(weights.ph2.unmod*B*(ii!=jj)*A.ij)), .as.double(V2), 
                    .as.double(dat.phase2[[Xcol]][jj]), as.integer(N), as.integer(dat[[dcol]]), .as.double(dat[[Xcol]]))
                V2 = -4 * matrix(tmp[[4]], nrow=p)
                
            }
            
            inv.gr = solve(gr)
            out = inv.gr %*% tcrossprod(V1+V2, inv.gr)
            
            # assuming that init is the truth, evaluate statistics at estimated survival function and true theta
            if (intermediate) attr(out, "intermediate") = cbind( calc(init, 4)/N**1.5, sqrt(diag(V1/N**3)), sqrt(diag((V1+V2)/N**3)) )
            # if we evaluate statistics at estimated survival function and estimated theta, then it is 0
            #if (intermediate) attr(out, "intermediate") = cbind( calc(theta, 4)/N**1.5, sqrt(diag(V1/N**3)), sqrt(diag((V1+V2)/N**3)) )
            #if (intermediate) attr(out, "intermediate") = cbind( colSums(h.rho.minus.1)/N**1.5, sqrt(diag(V1/N**3)), sqrt(diag((V1+V2)/N**3)) )#check h.rho.minus.1
            
        } else if (type==3) {
            #### compute efficient influence function
            out = -U %*% t(solve(gr))
            
        } else if (type==4) {
            #### compute a statistics useful for intermediate values
            out = colSums(U)
            
        } else if (type==5) {
            #### Kong et al 2013/8
            V1=matrix(0,p,p)
            tmp = .C("V1", as.integer(n), as.integer(p), .as.double((U)), .as.double(V1))
            V1 = matrix(tmp[[4]], nrow=p)
            
            
        } else {
            stop("type not supported")
        }
        
        return (out) 
    }
    
    
    # choosing abstol between 1, .1 and .01 does seem to make a big difference, abstol is used more than once
    abstol=0.01 
    # optim option trace=0 does nothing to the output
    optim.out=try(optim(init, calc, control=list(abstol=abstol, reltol=1e-100, maxit=maxit, trace=0)))
    
    est = rep(NA, p)
    if (!inherits(optim.out, "try-error")) {
        # when it does not converge, return NA
        if (optim.out$value<abstol) {
            est = optim.out$par 
        } 
    }
        
    #### get variance when est is not NA and we don't need efficient influence function
    #### when the intermediate attribute of .var is not null, also return it
    if(!any(is.na(est)) & !return.eif) {
        # assuming that init is the truth, evaluate statistics at estimated survival function and true theta
        #.var = calc(init, 2)
        .var = calc(est, 2)
        .sd = sqrt(diag(.var))
        tmp=attr(.var, "intermediate")
    } else {
        # if using NULL below, cbind with est will not create a second column
        .sd = NA
        tmp=NULL
    }
    out = cbind(est,"sd"=.sd)    
    if (!is.null(tmp)) out=cbind(out, tmp)
        
    #### this block attaches an attribute to out, so we shouldn't make certain changes to out after this
    #### evaluate efficient influence function at theta
    if(return.eif & !any(is.na(est))) {
        attr(out, "eif")=calc(est, 3)
    } else {
        # if using NA below, checking for NA when return.eif is present leads to a warning that there are more than one elements to compare
        attr(out, "eif")=NULL
    }
    
    # this may not be enough to release all memory
    gc()
    
    end=Sys.time()
    if (verbose) cat ("stm coef: ", out[,1], "\n")
    if(show.time.elapsed) cat("Time spent in stm.internal() "%+%format(end-begin)%+%"\n")
        
    attr(out,"class")="stm"    
    invisible (out) 
    
}


getFixedEf.stm=function (object, ...) {
    object[-1,1:2,drop=FALSE]
}
