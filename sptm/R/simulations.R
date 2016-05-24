# family="PH"; censoring.setting="1"; design="CC"; auxiliary="weak"; seed=1; var.S=1; var.W=1
# baseline.hazard 0.032
sim.fong = function (n, family=c("PH","PO","P2"), beta, 
    random.censoring=c("0%","20%","60%"), prevalence=0.1, non.adherence.ratio=0, 
    design=c("FULL","CC"), auxiliary=c("weak","good","excellent","none"), 
    seed=NULL, var.S=1, var.W=1) {
    
    random.censoring <- match.arg(random.censoring)
    design <- match.arg(design)    
    auxiliary <- match.arg(auxiliary)    
    
    if(!is.null(seed)) set.seed(seed)
    
    if (random.censoring=="0%") {
        C=rep(3,n) 
    } else if (random.censoring=="20%") {
        C=pmin(3, runif(n, 0, 15))
    } else if (random.censoring=="60%") {
        C=pmin(3, runif(n, 0, 5))
    }
    baseline.hazard = (prevalence-0.01133)/2.70362 # this is based on empirical results
    rho=switch(EXPR=auxiliary, weak=sqrt(.2), good=sqrt(.5), excellent=sqrt(.8), none=0)
        
    # covariates: z, s, w
    z = rep(c(0,1), each=n/2)
    sigma=diag(c(var.S, var.W))
    sigma[1,2] <- sigma[2,1] <- rho * sqrt(var.S * var.W)
    tmp = mvtnorm::rmvnorm(n, mean=c(0,0), sigma)
    s=tmp[,1]; w=tmp[,2]
    
    #### simulate failure time    
    #design.mat=model.matrix(formula, as.data.frame(cbind(X=1,d=1,z,s)))[,-1] # need dummy X, d; [,-1] removes intercept # not sure why it is done this way, replace with the following line
    design.mat=cbind(z,s,"z:s"=z*s)
    ft=rstm (n, family, design.mat %*% beta, baseline.hazard) 
    X=pmin(ft, C)
    d=ifelse(C>ft, 1, 0) # d = 1 means an event
    
    #### additional operations on S
    # loss due to non-adherence
    non.adherence.ind = rbinom(n, 1, non.adherence.ratio)
    s[non.adherence.ind==1]=NA
#    # early infection and early dropout: subjects not at risk anymore at 6 month cannot have s measured
#    s=ifelse (X<0.5, NA, s)    

    
    if (design=="FULL") {
        # full cohort
        out= (data.frame(ft, C, X, d, z, s))
        
    } else if (design=="CC") {    
        #### case-cohort design with finite population sampling
        
        phase2.p=0.22        
        subcohort=sample(n, n*phase2.p) # finite population sampling, as in Kong et al. (2004) and the usual practice
        subcohort.ind = rep(0,n); subcohort.ind[subcohort]=1        
    
        phase2.ind = d | subcohort.ind
    
        # phase 2 subjects do not have S measured
        s[phase2.ind==0] = NA        
    
        out = data.frame(ft, C, X, d, z, s, w) 
        
    } 
    
    invisible (out)
}


rstm=function (n, family=c("PH","PO","P2"), linear.predictors, baseline.hazard=1) {
    
    family <- match.arg(family)    
    exp.linear.predictors = exp(drop(linear.predictors)) # need to drop b/c linear.predictors is often obtained from matrix product
    x=runif(n)
    if (family=="PH") {
        ft = -log(1-x) / (exp.linear.predictors*baseline.hazard)
    } else if (family=="PO") {
        ft = log( 1+x/(1-x)/exp.linear.predictors ) / baseline.hazard
    } else if (family=="P2") {
        ft = log( 1+ (1/(1-x)^2-1)/exp.linear.predictors ) /2/baseline.hazard
    } 
    ft
}


sim.kong = function (gamma, beta, design="FULL", rho=0.9, seed=1, impute=FALSE, ppi) {
    
    set.seed(seed)
    
    n=1e3
    
    # simulate covariate
    x1=runif(n)
    x2=rbinom(n,1,.5)
    
    # simulate failure time
    ft=rstm (n, family, cbind(x1,x2) %*% beta, 1) 
    
    # simulate censoring time
    ub.c=.185 # for 90% censoring
    C=runif(n, 0, ub.c)
    d=ifelse(C>ft, 1, 0)
    X=pmin(ft, C)
    print("censoring proportion: " %+% round(1-mean(d),2))
    
    out=data.frame(ft, C, X, d, x1, x2) 
    
    if (design=="FULL") {
        indicators=rep(1,n)
    } else {
        # case-cohort design
        if (design=="CCI") subcohort.p = .11 else { if (design=="CCII") subcohort.p = .22 else stop("design not supported") }
        attr(out, "phase2.p")=subcohort.p
                
        out$indicators=d
        out$indicators[sample(n, n*subcohort.p)]=1
        print("control to case ratio " %+% round(sum(out$indicators)/sum(d)-1,2))
        
        # generate w
        if (rho==1) {
            out$w = logit(x1) 
        } else if (rho=="excellent") { # rho is approximately 0.9, not really sure 
            out$w = rnorm(n=n,sd=sqrt(1.72))+logit(x1) 
        } else if (rho=="good") {
            out$w = rnorm(n=n,sd=sqrt(7.2))+logit(x1) 
        } else if (rho==0 | rho=="NA") {
            out$w = rnorm(n=n,sd=sqrt(17.2))
        } else {
            stop ("are you sure rho is correct? ")
        }
#        mypostscript(mfrow=c(1,1))
#            plot(out$w[out$indicators==0], logit(out$x1)[out$indicators==0])
        
        if (!impute) {
            out[out$indicators==0,"x1"]=NA
            
            #### additional operation on x1
#            # early infection and early dropout: subjects not at risk anymore at 6 month cannot have s measured
#            x1=ifelse (X<0.05, NA, x1)    
#            print ("x1 is.na proportion: "%+%mean(is.na(x1)))
            # 15% loss due to non-adherence
            nonadhere=rbinom(nrow(out), 1, .15)==1            
            adherence.ratio = sum(out$indicators & !nonadhere) / sum(out$indicators)
            attr(out, "adherence.ratio")=adherence.ratio
            out[nonadhere,"x1"]=NA
            out[nonadhere,"indicators"]=0 # indicator should also reflect loss of s due to non-adherence and early dropout/infection
        } else {
            # impute x1 for all, not just those with missing values
            impute.fit = lm(logit(x1) ~ w, data=out, weights=1/ppi)
            out[,"x1"] = expit( predict(impute.fit,newdata=out[,"w",drop=F]) )
        }
    
#            lines(out$w[out$indicators==0], logit(out$x1)[out$indicators==0], col=2)
#        dev.off()
    }
    
    invisible( out )
}

# Simulate data as from Qin et al. Note that it may not work as baseline.hazard is not a global variable anymore. This function is currently not exported and called.
# rho.id=1; b3.id=1; missingness.id=1; BIP=TRUE; CPV=FALSE; seed=1
sim.qin = function (n, rho.id=1, b3.id=1, missingness.id=1, BIP=TRUE, CPV=FALSE, seed=1, var.S, var.W, b1, b2, b3, hazard, missingness, baseline.hazard) {
    
    set.seed(seed)
    
    b4=ifelse(BIP, 0, b4)
    rho=c(.9, .6)
    
    # simulate z
    z = rep(c(1,0), each=n/2)
    
    # simulate s and w
    sigma=diag(c(var.S, var.W))
    sigma[1,2]=sigma[2,1]= rho[rho.id] * sqrt(var.S * var.W)
    tmp = mvtnorm::rmvnorm(n, mean=c(0,0), sigma)
    s=tmp[,1]; w=tmp[,2]
    
    # compute lambda
    beta=matrix(c(b1, b2, b3[b3.id], b4), nrow=4, ncol=1)
    la = baseline.hazard * exp(cbind(z,s,z*s,w) %*% beta)
    
    # simulate failure time
    ft = drop( -log(1-runif(n))/la ) # failure time
    print(summary(ft))
    
    #    hist(ft[z==1], col=2, main="Failure time historgram", xlab="", xlim=c(0,200))
    #    hist(ft[z==0], add=T, col=3)
    #    legend(legend=c("vaccine","placebo"), col=2:3, lty=1, x="topright")
    #    abline(v=3)
    
    print(sum(ft[z==0]<3))
    print("should be around 334")
    
    # end of study censoring
    # d being 1 means an event happens
    d=rep(1,n)
    d[ft>3]=0
    X=ft
    X[X>3]=3
    
    # random censoring from 0 to a later boundary, so that mean ((ft>C)[ft<3]) is roughly 10% 
    C=runif(n, 0, 15)
#    d[ft>C]=0
#    # create X 
#    X=pmin(X,C)
    
    VAC = z==1
    UIV = z==1 & d==0
    IV  = z==1 & d==1 
    PLA = z==0
    UIP = z==0 & d==0
    IP  = z==0 & d==1 
    
    # remove some s that should be missing, weights are induced by missingness
    s.full=s
    probs=rep(1,n)
    
    strt=rep(1,n); strt[IV]=1; strt[UIV]=2; 
    
    # part of UIV missing due to two-phase sampling
    if (missingness[missingness.id]>0) s [which(UIV) [1:(sum(UIV)*missingness[missingness.id])]]=NA
    probs[UIV]=1-missingness[missingness.id]
    
    # UIP strata partially observed under CPV design
    if (CPV) {
        
        strt[IP]=3; strt[UIP]=4
        
        if (missingness[missingness.id]>0) s[which(UIP)[1:(sum(UIP)*missingness[missingness.id])]]=NA
        probs[UIP]=1-missingness[missingness.id]; 
    
        s.imputed.gold = s # at this stage, we have what we want to impute
    
        probs[IP]=0
        s[IP] = NA# IP strata should all be missing
        
        selected = !is.na(s)
        
        probs.imputed=probs; probs.imputed[IP]=1
        selected.imputed=selected; selected.imputed[IP]=TRUE
    } else {
        
        strt[PLA]=3
        
        s.imputed.gold = s # at this stage, we have what we want to impute
    
        s[PLA]=NA
        probs[PLA]=0
    
        selected = !is.na(s)
        
        probs.imputed=probs; probs.imputed[PLA]=1
        selected.imputed=selected; selected.imputed[PLA]=TRUE
    }
    
    dat.full         = data.frame(ft, C, X, d, z, w, strt, s=s.full)
    # the gold standard for imputation dataset
    dat.observed     = data.frame(ft, C, X, d, z, w, strt, s=s,              selected=selected,         probs=probs) 
    dat.imputed.gold = data.frame(ft, C, X, d, z, w, strt, s=s.imputed.gold, selected=selected.imputed, probs=probs.imputed) 
    #rm(ft, d, z, s, s.imputed.gold, w, strt, selected, selected.imputed, probs, probs.imputed, s.full)
    invisible(list(observed=dat.observed, full=dat.full, imputed.gold=dat.imputed.gold))
}


# ub.c is censoring distribution upper bound: 4 for 20-30% censoring, .4 for 80% censoring
sim.fine = function (n, ub.c, seed=1) {
    
    set.seed(seed)
    
    # simulate z
    x1=runif(n)
    x2=rbinom(n,1,.5)
    z = cbind(x1,x2)
    beta=matrix(c(-1,1), ncol=1)
    la = drop(exp(z %*% beta))
    
    # simulate failure time
    ft = -log(1-runif(n))/la  # failure time
    #print(summary(ft))
    #plot(density(ft))
    
    # simulate censoring time
    C=runif(n, 0, ub.c)
    d=C>ft
    X=pmin(ft, C)
    print("censoring proportion: " %+% (1-mean(d)))
        
    invisible( data.frame(ft, C, X, d, z) )
}


sim.kong1 = function (n, gamma, design="FULL", rho=0.9, seed=1, impute=FALSE, ppi) {
    
    set.seed(seed)
    # simulate z
    tmp=mvtnorm::rmvnorm(n, mean=c(0,0), sigma=matrix(c(0.4, 0.4*rho, 0.4*rho, 0.4),2))
    x1=tmp[,1]
    w=tmp[,2]
    x2=rbinom(n,1,.5)
    z = cbind(x1,x2)
    beta=matrix(c(1,-1), ncol=1)
    la = drop(exp(z %*% beta))
#    mypostscript(mfrow=c(1,2))
        plot(w, x1, xlim=c(-2,2), ylim=c(-2,2))
        abline(lm(x1 ~ w))
        summary(lm(x1 ~ w))
    
    # simulate failure time
    if (gamma==0) {
        ft = -log(1-runif(n))/la  # failure time
    }
    #print(summary(ft))
    #plot(density(ft))
    
    # simulate censoring time
    ub.c=.185 # for 90% censoring
    C=runif(n, 0, ub.c)
    d=ifelse(C>ft, 1, 0)
    X=pmin(ft, C)
    print("censoring proportion: " %+% round(1-mean(d),2))
    
    out=data.frame(ft, C, X, d, z, w) 
    
    if (design=="FULL") {
        indicators=rep(1,n)
    } else {
        # case-cohort design
        if (design=="CCI") subcohort.p = .11 else { if (design=="CCII") subcohort.p = .22 else stop("design not supported") }
        attr(out, "subcohort")=n*subcohort.p
        attr(out, "cohort")=n
        attr(out, "controls")=n-sum(d)
                
        out$indicators=d
        out$indicators[sample(n, n*subcohort.p)]=1
        attr(out, "selectedcontrols")=sum(out$indicators)-sum(d)
        print("control to case ratio " %+% round(sum(out$indicators)/sum(d)-1,2))
        out$ppi = 1
        out$ppi[d==0] = subcohort.p
        #out$ppi[d==0] = attr(out, "selectedcontrols") / attr(out, "controls")
            
        if (!impute) {
            # only return complete cases
            out=out[out$indicators==1,]
        } else {
            # impute x1 for all, not just those with missing values
            actual=out$x1[out$indicators==1]
            impute.fit = lm(x1 ~ w, data=out, weights=1/ppi)
            print(summary(impute.fit))
            out[,"x1"] = predict(impute.fit,newdata=out[,"w",drop=F]) 
            imputed=out$x1[out$indicators==1]
    
            plot(actual, imputed, xlim=c(-2,2), ylim=c(-2,2))
            lines(lowess(imputed~actual))
            dev.off()
        }
    
    }
        
    invisible( out )
}
