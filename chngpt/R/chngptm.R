#deviance.3pl <- expression( (1-y) * ( c + (d-c)/(1+exp(b*(t-e))) )  +  log( 1 + exp( -c - (d-c)/(1+exp(b*(t-e))) ) ) )
#deviance.3pl.deriv=deriv3(deviance.3pl, c("c","d","e"), c("c","d","e","t","y","b"))

# step change point model
dev.step <- expression( (1-y) * ( alpha.z + beta/(1+exp(b*(x-e))) )  +  log( 1 + exp( -alpha.z - beta/(1+exp(b*(x-e))) ) ) )
dev.step.deriv=deriv3(dev.step, c("beta","e"), c("beta","e","x","y","b","alpha.z"))
dev.step.itxn <- expression( (1-y) * ( alpha.z + (beta1+beta2*z.1)/(1+exp(b*(x-e))) )  +  log( 1 + exp( -alpha.z - (beta1+beta2*z.1)/(1+exp(b*(x-e))) ) ) )
dev.step.itxn.deriv=deriv3(dev.step.itxn, c("beta1","beta2","e"), c("beta1","beta2","e","x","y","b","alpha.z","z.1"))
# hinge change point model
dev.hinge <- expression( (1-y) * ( alpha.z + (x-e)*beta/(1+exp(b*(x-e))) )  +  log( 1 + exp( -alpha.z - (x-e)*beta/(1+exp(b*(x-e))) ) ) )
dev.hinge.deriv=deriv3(dev.hinge, c("beta","e"), c("beta","e","x","y","b","alpha.z"))
dev.hinge.itxn <- expression( (1-y) * ( alpha.z + (x-e)*(beta1+beta2*z.1)/(1+exp(b*(x-e))) )  +  log( 1 + exp( -alpha.z - (x-e)*(beta1+beta2*z.1)/(1+exp(b*(x-e))) ) ) )
dev.hinge.itxn.deriv=deriv3(dev.hinge.itxn, c("beta1","beta2","e"), c("beta1","beta2","e","x","y","b","alpha.z","z.1"))

# step change point model
dev.step.f <- function(theta,x,y,b,alpha.z,z.1)  {
    beta=theta[1]; e=theta[2]
    eta=beta/(1+exp(b*(x-e))) 
    linear=alpha.z+eta
    sum( (1-y)*linear + log(1+exp(-linear)) )
}
dev.step.itxn.f <- function(theta,x,y,b,alpha.z,z.1)  {
    beta1=theta[1]; beta2=theta[2]; e=theta[3]
    eta=(beta1+beta2*z.1)/(1+exp(b*(x-e))) 
    linear=alpha.z+eta
    sum( (1-y)*linear + log(1+exp(-linear)) )
}
# hinge change point model
dev.hinge.f <- function(theta,x,y,b,alpha.z,z.1) {
    beta=theta[1]; e=theta[2]
    eta=(x-e)*beta/(1+exp(b*(x-e)))
    linear=alpha.z+eta
    sum( (1-y)*linear + log(1+exp(-linear)) )
} 
dev.hinge.itxn.f <- function(theta,x,y,b,alpha.z,z.1)  {
    beta1=theta[1]; beta2=theta[2]; e=theta[3]
    eta=(x-e)*(beta1+beta2*z.1)/(1+exp(b*(x-e))) 
    linear=alpha.z+eta
    sum( (1-y)*linear + log(1+exp(-linear)) )
}
# segmented change point model
dev.segmented.f <- function(theta,x,y,b,alpha.z,z.1) {
    beta1=theta[1]; beta2=theta[2]; e=theta[3]
    eta=beta1*x + (x-e)*beta2/(1+exp(b*(x-e))) 
    linear=alpha.z+eta
    sum( (1-y)*linear + log(1+exp(-linear)) )
} 
dev.segmented.itxn.f <- function(theta,x,y,b,alpha.z,z.1) {
    beta1=theta[1]; beta2=theta[3]; beta3=theta[2]; beta4=theta[4]; e=theta[5]# note the order change between beta and theta
    eta=(beta1+beta2*z.1)*x + (beta3+beta4*z.1)*(x-e)/(1+exp(b*(x-e)))
    linear=alpha.z+eta
    sum( (1-y)*linear + log(1+exp(-linear)) )
} 
# stegmented change point model
dev.stegmented.f <- function(theta,x,y,b,alpha.z,z.1) {
    beta1=theta[1]; beta2=theta[2]; beta3=theta[3]; e=theta[4]
    eta=beta1*x + ((x-e)*beta3 + beta2)/(1+exp(b*(x-e))) 
    linear=alpha.z+eta
    sum( (1-y)*linear + log(1+exp(-linear)) )
} 
# the order of parameters in the following needs to be fixed
dev.stegmented.itxn.f <- function(theta,x,y,b,alpha.z,z.1) {
    beta1=theta[1]; beta2=theta[3]; beta3=theta[2]; beta4=theta[4]; beta5=theta[5]; beta6=theta[6]; e=theta[5]# note the order change between beta and theta
    eta=(beta1+beta2*z.1)*x + ((beta3+beta4*z.1)*(x-e) + beta4+beta5*z.1)/(1+exp(b*(x-e)))
    linear=alpha.z+eta
    sum( (1-y)*linear + log(1+exp(-linear)) )
} 


do.regression=function(formula.new, data, weights, family){
    if (family=="coxph") {
        fit = keepWarnings(survival::coxph(formula.new, data=data, weights=weights))
    } else {
        fit =             keepWarnings(glm(formula.new, data=data, weights=weights, family=family) )
    }
    fit
}


# tol=1e-4; maxit=1e2; verbose=TRUE; est.method="grid"; lb.quantile=.1; ub.quantile=.9; grid.size=500
chngptm = function(formula.1, formula.2, data, family,
  type=c("step","hinge","segmented","stegmented"),
  est.method=NULL, 
  lb.quantile=.1, ub.quantile=.9, grid.size=500,
  tol=1e-4, maxit=1e2, chngpt.init=NULL, 
  weights=NULL,
  verbose=FALSE,
  ...) 
{
    
    # remove missing observations
    form.all = update(formula.1, formula.2)
    subset. = complete.cases(model.frame(form.all, data, na.action=na.pass))
    data=data[subset.,,drop=FALSE]
    
    if (missing(type)) stop("type mssing")
    type<-match.arg(type)    
    p.2=switch(type, step=1, hinge=1, segmented=2, stegmented=3)
    
    if (is.null(est.method)) 
        est.method = if(family=="binomial" & is.null(weights) & type=="step") "sigmoidapprox" else "grid"
    
    b.=-30    
        
    # decide whether there is interaction
    y=model.frame(formula.1, data)[,1]
    Z=model.matrix(formula.1, data)
    tmp=model.matrix(formula.2, data)[,-1,drop=F]
    n=nrow(Z)
    p=ncol(Z)
    chngpt.var.name=setdiff(colnames(tmp), colnames(Z))[1]
    z.1.name=intersect(colnames(tmp), colnames(Z))
    chngpt.var = tmp[,chngpt.var.name]
    z.1 = tmp[,z.1.name] # if the intersection is a null set, z.1 is a matrix of n x 0 dimension
    has.itxn = length(z.1.name)>0
    p.2.itxn=p.2*ifelse(has.itxn,2,1)
    
    if (is.null(weights)) weights=rep(1,n) 
    data$weights=weights # try put it in data to be found by glm
    
    # make formula
    if (type %in% c("segmented","stegmented")) formula.new = update(formula.1, as.formula("~.+"%+%chngpt.var.name)) else formula.new=formula.1
    f.alt=get.f.alt(type, has.itxn, z.1.name, chngpt.var.name)
    formula.new=update(formula.new, as.formula(f.alt))
    
    if (verbose) myprint(type, est.method, has.itxn)     
    if (verbose) print(formula.new)
        
    if (est.method=="grid") {
        
        # define a list of potential change points
        chngpts=sort(chngpt.var)
        chngpts=chngpts[chngpts<quantile(chngpts, ub.quantile) & chngpts>quantile(chngpts, lb.quantile)]
        if (length(chngpts)>grid.size) chngpts=chngpts[round(seq(1,length(chngpts),length=grid.size))]
                
        logliks=sapply (chngpts, function(e) {
            data=make.chngpt.var(chngpt.var, e, type, data)
            fit = do.regression (formula.new, data, weights, family)
            if(length(fit$warning)!=0) {
                if(verbose) print(fit$warning)
                NA
            } else as.numeric(logLik(fit$value))
        } )
        glm.warn=any(is.na(logliks))
        e=chngpts[which.max(logliks)]
        if(verbose==2) {
            plot(chngpts, logliks, type="b", xlab="change points")
            abline(v=e)
        }
        data = make.chngpt.var(chngpt.var, e, type, data)
        fit = do.regression (formula.new, data, weights, family)$value
        coef.hat=c(coef(fit), "chngpt"=e)       
        
    } else if (est.method=="sigmoidapprox") {
    
        # Newton-Raphson
        
        # do a test to get init value for change point to be used in estimation
        if (is.null(chngpt.init)) {
            test = chngpt.test (formula.1, formula.2, data, type=type) # there is no test for segmented or itxn, yet
            if (verbose==2) print(test)
            if (verbose==2) plot(test)
            e.init=test$chngpt
        } else {
            test=NULL # test is returned either way
            e.init=chngpt.init
        }
        names(e.init)="e"
        glm.warn=test$glm.warn
        
        coef.hat=rep(0, 1+ncol(Z)+p.2.itxn)
        n.iter=0
        converged=TRUE        
        while(TRUE){    
        
            n.iter=n.iter+1
            if (verbose) cat("iter ", n.iter, "\n")
            if (n.iter>maxit) {converged=FALSE; break}
            
            # remake the binary change point variable in every iteration based on the change point estimate from last iteration
            data = make.chngpt.var(chngpt.var, e.init, type, data)            
            fit.0 = do.regression(formula.new, data, weights, family)$value
            if (verbose==2) print(coef(fit.0))
            
            # update threshold and associated coefficients
            beta.init=coef(fit.0)[p+1:p.2.itxn]; #names(beta.init)=c("beta1","beta2")# we may need this for variance calculation to work
            alpha.hat=coef(fit.0)[1:p]
            alpha.z = c(Z %*% alpha.hat)
            if (verbose) myprint(beta.init, e.init, b.)
    
            optim.out = optim(par=c(beta.init, e.init), 
                  fn = get("dev."%+%type%+%"."%+%ifelse(has.itxn,"itxn.","")%+%"f"), 
                  # if we use analytical gradient function by deriv3 in optim, we can get situations like exp(100), which will be Inf, and Inf/Inf will be NaN
                  gr = NULL,
                  #fn = function(theta,...) sum(dev.step.itxn.deriv(theta[1],theta[2],theta[3],...)), 
                  #gr = function(theta,...) colSums(attr(dev.step.itxn.deriv(theta[1],theta[2],theta[3],...), "gradient")), 
                  chngpt.var, y, b., alpha.z, z.1,
                  lower = c(rep(-10,length(beta.init)), quantile(chngpt.var, lb.quantile)), 
                  upper = c(rep( 10,length(beta.init)), quantile(chngpt.var, ub.quantile)), 
                  method="L-BFGS-B", control = list(), hessian = TRUE)
            
            e.init=optim.out$par["e"]; names(e.init)="e"
            coef.tmp=c(alpha.hat, optim.out$par)
            if (verbose) print(coef.tmp)
            if (max(abs(coef.tmp-coef.hat))<tol) {
                coef.hat=coef.tmp
                break
            } else {
                coef.hat=coef.tmp
            }
            
        } # end while 
    
    } # end if est.method 
    
    if (verbose) print(coef.hat)
    names(coef.hat)[length((coef.hat))]="chngpt"
    if (type=="stegmented") {
        replacement="I("%+%chngpt.var.name%+%">chngpt)"
        replacement.2="I("%+%chngpt.var.name%+%">chngpt)*("%+%chngpt.var.name%+%"-chngpt)"
        new.names=sub("x.gt.e.2", replacement.2, names(coef.hat))    
        new.names=sub("x.gt.e", replacement, new.names)    
    } else {
        if(type=="step") {
            replacement="I("%+%chngpt.var.name%+%">chngpt)"
        } else if (type=="hinge" | type=="segmented") {
            replacement="I("%+%chngpt.var.name%+%">chngpt)*("%+%chngpt.var.name%+%"-chngpt)"
        } 
        new.names=sub("x.gt.e", replacement, names(coef.hat))    
    }
    if (verbose) print(new.names)
    names(coef.hat)=new.names
    
#    # variance-covariance matrix
#    if (type=="step" & family=="binomial") {
#        # expressions for use with optim with deriv3
#        alpha.z.s="("%+% concatList("alpha"%+%1:p%+%"*z"%+%1:p%+%"","+") %+%")"
#        if (!has.itxn) {
#            deviance.s <- " (1-y) * ( "%+% alpha.z.s %+%" + beta/(1+exp(b*(x-e))) )  +  log( 1 + exp( -"%+% alpha.z.s %+%" - beta/(1+exp(b*(x-e))) ) ) "
#            params=c("alpha"%+%1:p, "beta", "e")
#        } else {
#            deviance.s <- " (1-y) * ( "%+% alpha.z.s %+%" + (beta1+beta2*z.1)/(1+exp(b*(x-e))) )  +  log( 1 + exp( -"%+% alpha.z.s %+%" - (beta1+beta2*z.1)/(1+exp(b*(x-e))) ) ) "
#            params=c("alpha"%+%1:p, "beta1", "beta2", "e")
#        }
#        params.long=c(params,"x","y","b","z"%+%1:p,"z.1")
#        if (verbose) print(deviance.s)
#        if (verbose) myprint(params)
#        loss.f=deriv3(parse(text=deviance.s), params, params.long)    
#        param.list = c(as.list(coef.hat), list(chngpt.var), list(y), list(b.), lapply(1:ncol(Z), function (i) Z[,i]), list(z.1))
#        names(param.list)=params.long    
#        tmp=do.call(loss.f, param.list)
#        hess=apply(attr(tmp,"h"), 2:3, sum, na.rm=T)                
#        var.est = try(solve(hess)) # should keep change point in, and not do hess[-ncol(hess), -ncol(hess)], otherwise lead to over estimation of sd
#        rownames(var.est) <- colnames(var.est) <- names(coef.hat)
#    } else {
#        var.est=NULL
#    }
    
    res=list(
          coefficients=coef.hat
#        , vcov=var.est
        , formula.1=formula.1
        , formula.2=formula.2
        , chngpt.var=chngpt.var.name
        , chngpt=coef.hat["chngpt"]
        , est.method = est.method
        , glm.warn = glm.warn    
    )
    if (est.method=="sigmoidapprox") {
        res=c(res, list(converged=converged, iter=n.iter, test=test, best.fit=fit.0))
    } else if (est.method=="grid") {
        res=c(res, list(best.fit=fit))
    }
    names(res$chngpt) = round(100*mean(chngpt.var<res$chngpt),1) %+% "%"
    class(res)=c("chngptm", class(res))
    
    res    
}

print.chngptm=function(x, ...) {
    if (x$est.method=="sigmoidapprox") {
        if (!x$converged) cat("Warning: not converged\n")
    }
    print(x$coefficients)
}

coef.chngptm=function(object, ...) {
    object$coefficients[-length(object$coefficients)] # not return the chngpoint estimate
}

vcov.chngptm=function(object, ...) {
    object$vcov[-length(object$coefficients),-length(object$coefficients)] # not return the chngpoint estimate
}

summary.chngptm=function(object, ...) {
    
    # "Estimate" "Std. Error" "t value" "Pr(>|t|)"    
    
    fit=object
    p=length(fit$coefficients)
    n=nrow(fit$data)
    
    if (is.null(fit$vcov)) {
        cat("No variance estimate available.\n\n")
        print(fit)
        return (invisible())
    }
    
    # assuming the last of coefficients is always the change point
    res=list()
    res$coefficients=mysapply(1:(p-1), function (i) {
        c(
              "OR"=exp(unname(fit$coefficients[i]))
            , "p.value" = unname(pt(abs(fit$coefficients[i] / sqrt(fit$vcov[i,i])), df=n-p, lower.tail=FALSE))
            , "(lower" = exp(unname(fit$coefficients[i] - sqrt(fit$vcov[i,i]) * qt(0.975, df=n-p, lower.tail=TRUE)))
            , "upper)" = exp(unname(fit$coefficients[i] + sqrt(fit$vcov[i,i]) * qt(0.975, df=n-p, lower.tail=TRUE)))
#              "Estimate"=unname(fit$coefficients[i])
#            , "Std. Error" = sqrt(fit$vcov[i,i])
#            , "t value" = unname(fit$coefficients[i] / sqrt(fit$vcov[i,i]))
#            , "Pr(>|t|)" = unname(pt(abs(fit$coefficients[i] / sqrt(fit$vcov[i,i])), df=n-p, lower.tail=FALSE))
        )
    })
    rownames(res$coefficients)=names(fit$coefficients)[-p]
    
    i=p
    res$chngpt=c(
              fit$chngpt
            , "p.value" = unname(pt(abs(fit$coefficients[i] / sqrt(fit$vcov[i,i])), df=n-p, lower.tail=FALSE))
            , "(lower" = (unname(fit$coefficients[i] - sqrt(fit$vcov[i,i]) * qt(0.975, df=n-p, lower.tail=TRUE)))
            , "upper)" = (unname(fit$coefficients[i] + sqrt(fit$vcov[i,i]) * qt(0.975, df=n-p, lower.tail=TRUE)))
    )
    
    res
}


#    c.init=ifelse(beta<0, (alpha+beta), (alpha)); names(c.init)="c"
#    d.init=ifelse(beta>0, (alpha+beta), (alpha)); names(d.init)="d"            
#    lb=c(logit(1/n),logit(1/n),quantile(dat$x,.01))
#    ub=c(logit((n-1)/n),logit((n-1)/n),quantile(dat$x,.99))
#        
#    ee=quantile(dat$x,seq(.15,.85,length=10)); ee=c(ee, e.)
#    fits=lapply(ee, function (e.init){
#        names(e.init)="e"
#        optim(par=c(c.init,d.init,e.init), 
#              fn = function(theta,...) sum(deviance.3pl.deriv(theta[1],theta[2],theta[3],...)), 
#              gr = function(theta,...) colSums(attr(deviance.3pl.deriv(theta[1],theta[2],theta[3],...), "gradient")), 
#              dat$x, dat$y, b., method="L-BFGS-B", lower = lb, upper = ub, control = list(), hessian = F)
#    })            
#    fit=fits[[which.min(sapply(fits, function(x) x$value))]]
#    res=c(res, "fit"=chngpt.score.stat(e.=fit$par["e"], y~x, dat, b.=b.)["Z.stat"]) # similar performance to max.Z
#    
#    
#    plot(y~x, dat, log=""); abline(v=exp(e.))
#    
#    # likelihood evaluated at init
#    tmp=sum(deviance.3pl.deriv(c.init,d.init,e.init,dat$x,dat$y)); sum(tmp)
#    
#    # plotting likelihood
#    ee=exp(seq(3,7,length=100))
#    llik=sapply(ee, function (e.1){
#        lik=ThreePL.chng.t(t=dat$x,c=c.init,d=d.init,e=e.1)
#        lik[dat$y==0]=1-lik[dat$y==0]
#        llik=sum(log(lik))
#    })
#    plot(ee, llik, log="x", type="l")
#    
#    lik=ThreePL.chng.t(t=dat$x,c=coef(fit)[1],d=coef(fit)[2],e=coef(fit)[3])
#    lik[dat$y==0]=1-lik[dat$y==0]
#    sum(log(lik))
    
#    # compare variance estimate computed by deriv3 and by optim
#    tmp=deviance.3pl.deriv(coef(fit)[1],coef(fit)[2],coef(fit)[3],dat$x,dat$y)
#    sum(tmp, na.rm=T)
#    solve(-apply(attr(tmp,"h"), 2:3, mean))
#    
#    fit$fit$val
#    solve(fit$fit$hessian) # same as vcov(fit)
#    
#    
#    if (fit$convergence!=0) { 
#        did.not.fit=T
#    } else {
#        # compute variance estimate using deriv3
#        tmp=deviance.3pl.deriv(fit$par["c"], fit$par["d"], fit$par["e"], dat$x, dat$y)
#        hess=apply(attr(tmp,"h"), 2:3, sum, na.rm=T)                
#        var.est = try(solve(hess)[1:2,1:2])
#        if (class(var.est)=="try-error") {
#            did.not.fit=T
#        } else {
#            did.not.fit=F
#        }
#    }
#    
#    if (did.not.fit) {
#        res=c(res, sigm=NA, sigm.cover=NA)
#    } else {
#        d.minus.c = fit$par["d"] - fit$par["c"]
#        d.minus.c.var = c(-1,1) %*% var.est %*% c(-1,1) 
#        res=c(res, sigm = pt(abs(d.minus.c)/sqrt(d.minus.c.var), df=nrow(dat)-3, lower.tail=F) * 2 )                                 
#        res=c(res, sigm.cover = d.minus.c - 1.96*sqrt(d.minus.c.var) < beta & beta < d.minus.c + 1.96*sqrt(d.minus.c.var) ) # coverage
#        res=c(res, fit$par["e"]) # e
#    }
