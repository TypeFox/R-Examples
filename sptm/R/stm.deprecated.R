###########################################
# deprecated
# p.weight: probability weight, should be "", "IPW", "BRESLOW", if it is a vector, it should have length n^2, where n is nrow(dat)
stm.deprecated=function(formula, dat, interest.model, method, init, p.weight="", t0=NULL, eff.inf=F) {
    
    # remove rows with missing data
    vars.in.formula=intersect(all.vars(interest.model), attr(terms(interest.model),"term.labels"))
    dat=dat[complete.cases(dat[,vars.in.formula]),] # using subset removes attributes , attr
    
    # create building blocks for estimating functions
    Z=model.matrix(formula, dat)[,-1,drop=F] # remove intercept column
    if (nrow(dat)!=nrow(Z)) stop ("nrow(dat)!=nrow(Z)")    
    n=nrow(Z)
    ii = rep (1:n, n)
    jj = rep (1:n, each=n)    
    Z.ij = Z[ii,,drop=F]-Z[jj,,drop=F] # n**2 by p
    Z.ii = Z[ii,,drop=F] # n**2 by p
    Z.jj = Z[jj,,drop=F] # n**2 by p
    
    weights=rep(1,n*n)
    if (length(p.weight)>1) {
        weights = p.weight
    } else if (p.weight=="") {
        # do nothing 
    } else if (p.weight=="IPW") {
        subcohort=attr(dat, "subcohort")
        cohort=attr(dat, "cohort")
        weights[dat$d[ii]==1 & dat$d[jj]==0]=subcohort/cohort
        weights[dat$d[ii]==0 & dat$d[jj]==1]=subcohort/cohort
        weights[dat$d[ii]==0 & dat$d[jj]==0]=subcohort/cohort*(subcohort-1)/(cohort-1)
        weights=1/weights
    } else if (p.weight=="IPW1") {
        subcohort=attr(dat, "subcohort")
        cohort=attr(dat, "cohort")
        weights[dat$d[ii]==1 & dat$d[jj]==0]=subcohort/cohort
        weights[dat$d[ii]==0 & dat$d[jj]==1]=subcohort/cohort
        weights[dat$d[ii]==0 & dat$d[jj]==0]=(subcohort/cohort)**2
        weights=1/weights
    } else if (p.weight=="EST") {
        est.p=attr(dat, "selectedcontrols")/attr(dat, "controls")
        weights[dat$d[ii]==1 & dat$d[jj]==0]=est.p
        weights[dat$d[ii]==0 & dat$d[jj]==1]=est.p
        weights[dat$d[ii]==0 & dat$d[jj]==0]=est.p*est.p
        weights=1/weights
    } 
    
    # KM estimate of the censoring distribution
    censoring.s = survfit(Surv(X,1-d)~1, data=dat)
    #plot(censoring.s)    
    # old
#    idx=sapply(1:n, function(j) sum(censoring.s$time<dat$X[j]) )
#    G.hat = ifelse(idx==0, 1, censoring.s$surv[idx])
    # new
    idx=sapply(1:n, function(j) sum(censoring.s$time<=dat$X[j]) )
    G.hat = censoring.s$surv[idx]
    #G.hat = G.hat + .00001 # to avoid the problem of dividing by 0
    
    if (startsWith(method, "fine") | method=="kong" ) {
        # truncation point
        if (is.null(t0)) t0=quantile(dat[dat$d==1,"X"],.85) # Fine et al choice, default
        cat("truncation point: " %+% t0, "\n    ")
        observed = with(dat, ifelse(t0>=X[jj], d[jj] * (X[ii]>=X[jj]) / G.hat[jj]**2, 0) ) # d and X are part of dat
    }
    if (method=="cheng") {
        observed = with(dat, d[jj] * (X[ii]>=X[jj]) / G.hat[jj]**2) # d and X are part of dat
    }
    
    ee.kong=function(param){
        alpha=param[1]; beta=matrix(param[-1], ncol=1)
                        
        e.ii = c( exp(alpha+Z.ii%*%beta) )
        e.jj = c( exp(alpha+Z.jj%*%beta) )
        e.ij = c( exp(Z.ij%*%beta) )
        e.nij = 1/e.ij
        
        expected = 1/(1+e.ij)*(1-exp(-e.ii-e.jj))
        obs.minus.exp = weights * (ii!=jj) * (observed-expected) / expected / (1-expected)
        
        a1 = 1/(1+e.ij)-expected            
        ee.alpha = sum ( obs.minus.exp * a1 * (e.ii + e.jj) )
        ee.beta = colSums ( obs.minus.exp * ( a1 * (Z.ii*e.ii+Z.jj*e.jj) - expected/(1+e.nij) * Z.ij) )
        out = ee.alpha**2 + sum(ee.beta**2)
        out
    }
    
    ee.fine=function(param){
        alpha=param[1]; beta=matrix(param[-1], ncol=1)
                        
        e.ii = c( exp(alpha+Z.ii%*%beta) )
        e.jj = c( exp(alpha+Z.jj%*%beta) )
        e.ij = c( exp(Z.ij%*%beta) )
        e.nij = 1/e.ij
        
        expected = 1/(1+e.ij)*(1-exp(-e.ii-e.jj))
        obs.minus.exp = weights * (ii!=jj) * (observed-expected)
        
        a1 = 1/(1+e.ij)-expected            
        ee.alpha = sum ( obs.minus.exp * a1 * (e.ii + e.jj) )
        ee.beta = colSums ( obs.minus.exp * ( a1 * (Z.ii*e.ii+Z.jj*e.jj) - expected/(1+e.nij) * Z.ij) )
        out = ee.alpha**2 + sum(ee.beta**2)
        out
    }
    
    # fine1 differs from fine in ee.beta, this is the fastest of fine variations
    ee.fine1=function(param){
        alpha=param[1]; beta=matrix(param[-1], ncol=1)
                        
#        e.ii = c( exp(alpha+Z.ii%*%beta) )
#        e.jj = c( exp(alpha+Z.jj%*%beta) )
#        e.ij = c( exp(Z.ij%*%beta) )
        
        a1 = 1+ c( exp(Z.ij%*%beta) )
        a2 = c( exp(alpha+Z.jj%*%beta) ) * a1
                
        expected = 1/a1*(1-exp(-a2))
        obs.minus.exp = weights * (ii!=jj) * (observed-expected)
        
        ee.alpha = sum ( obs.minus.exp * (1/a1-expected) * a2 )
        ee.beta = colSums ( obs.minus.exp * ( Z.ij ) )
        out = ee.alpha**2 + sum(ee.beta**2)
        
        print(out)
        out
    }
    
    # fine2 differs from fine in ee.beta and ee.alpha, this is the simplest of fine variations
    ee.fine2=function(param){
        alpha=param[1]; beta=matrix(param[-1], ncol=1)
                        
        e.ii = c( exp(alpha+Z.ii%*%beta) )
        e.jj = c( exp(alpha+Z.jj%*%beta) )
        e.ij = c( exp(Z.ij%*%beta) )
        
        expected = 1/(1+e.ij)*(1-exp(-e.ii-e.jj))
        obs.minus.exp = weights * (ii!=jj) * (observed-expected)
        
        ee.alpha = sum ( obs.minus.exp * 1 )
        ee.beta = colSums ( obs.minus.exp * ( Z.ij ) )
        out = ee.alpha**2 + sum(ee.beta**2)
        #cat(out," ")
        out
    }
    # gradient function
    gr.fine2=function (param){
        alpha=param[1]; beta=matrix(param[-1], ncol=1)
        # c() is necessary so that row-wise operation can be performed between, say, nx2 matrix and nx1 matrix 
        a1 = 1/(1 + c( exp(Z.ij%*%beta) ))
        a2 = c( exp(alpha+Z.jj%*%beta) ) / a1
        a3 = Z.ij / c(1+exp(-Z.jj%*%beta))
        gr=-t(cbind(1, Z.ij)) %*% ((ii!=jj)*a1*cbind( a2*exp(-a2), exp(-a2)*( Z.ii*c(exp(alpha+Z.ii%*%beta)) + Z.jj*c(exp(alpha+Z.jj%*%beta)) + a3 )-a3))
    }
    # efficient influence function
    ei.fine2=function(param) {
        # copied from gr.xxxx function, because we don't want to compute a1 and a2 twice
        alpha=param[1]; beta=matrix(param[-1], ncol=1)
        # c() is necessary so that row-wise operation can be performed between, say, nx2 matrix and nx1 matrix 
        a1 = 1/(1 + c( exp(Z.ij%*%beta) ))
        a2 = c( exp(alpha+Z.jj%*%beta) ) / a1
        a3 = Z.ij / c(1+exp(-Z.jj%*%beta))
        gr=-t(cbind(1, Z.ij)) %*% ((ii!=jj)*a1*cbind( a2*exp(-a2), exp(-a2)*( Z.ii*c(exp(alpha+Z.ii%*%beta)) + Z.jj*c(exp(alpha+Z.jj%*%beta)) + a3 ) -a3))
    
        expected = a1*(1-exp(-a2))
        U = (ii!=jj) * (observed-expected) * cbind(1,Z.ij)
        
        ei = -U %*% t(solve(gr))
    }
    
    ee.cheng=function(param){
        beta=matrix(param, ncol=1)
                        
        e.ij = c( exp(Z.ij%*%beta) )
        
        expected = 1/(1+e.ij)
        obs.minus.exp = weights * (ii!=jj) * (observed-expected)
        
        out=sum(colSums(obs.minus.exp * Z.ij)**2)        
        out
    }
    
    # init has to be assigned here because it relies on t0 which may have to be assigned in this function 
    #if (method!="cheng") init=c(log(t0*baseline.hazard),init)
    if (method!="cheng") init[1]=init[1]+log(t0)
    
    ee=get("ee."%+%method)
    if(eff.inf) ei=get("ei."%+%method)
    
    # using cox estimate as init does not help.
    # choosing this between 1, .1 and .01 does seem to make a big difference
    # optim option trace=0 does nothing to the output
    maxit=1000 # necessary for simulation "fong" because it may run for a long time otherwise
    abstol=0.01 # this variable is used more than once
    tmp=optim(init, ee, control=list(abstol=abstol, reltol=1e-100, maxit=maxit))
    print("") # to add a line break
    print (str(tmp))
    if (tmp$value<abstol) {
        est = tmp$par 
    } else {
        est = rep(NA, length(tmp$par))
#        print("change an init value to do optim again")
#        tmp=optim(init*3, ee, control=list(abstol=abstol, reltol=1e-100))
#        print(tmp)
#        if (tmp$value<abstol) {
#            est = tmp$par 
#        } else {
#            est = rep(NA, length(tmp$par))
#        }
    }
    print(est)
        
    if(eff.inf & !any(is.na(est))) {
        # get efficient influence function and compute raking weights
        r=ei(est)
        xi=dat$indicators[ii] & dat$indicators[jj]
        selection.prob = dat$ppi[ii] * dat$ppi[jj]
        raking=function(lambda) {
            est.func = c(xi / selection.prob * exp(- r %*% lambda) - weights ) * r
            out = sum( colSums(est.func) ** 2 )
            #cat(out," ")
            out
        }
        print("compute raking weights ...")
        tmp=optim(rep(0,ncol(r)), raking, control=list(abstol=abstol, reltol=1e-100, maxit=maxit))
        print("") # to add a line break
        print (str(tmp))
        if (tmp$value<abstol) {
            lambda = tmp$par 
        } else {
            print ("raking fail to converge, return NA ... ")
            return (rep(NA, length(tmp$par)))
        }
        print(lambda)
        raking.w = ( exp(- r %*% lambda) / selection.prob )[xi]
        return (raking.w)
    } else {
        return (est)
    }
        
} 
