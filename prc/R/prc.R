# reltol=1e-3; opt.method="gnls"; max.iter=50; verbose=TRUE; model="4P"; method="TLS"; init.method="gnls"
prc=function(xvar, dil.x, yvar, dil.y, model=c("4P","3P"), method=c("TLS","naive"),
  init.method=c("gnls","optim"), opt.method=c("gnls","optim"), reltol=1e-3, max.iter=50, init=NULL,
  verbose=FALSE) {    
    
    model=match.arg(model)
    method=match.arg(method)
    opt.method=match.arg(opt.method)
    init.method=match.arg(init.method)
    
    if (verbose) {myprint(model, method, opt.method, init.method)}
    
    f.is.1 = model=="3P"
    if(f.is.1 & opt.method!="gnls") stop("3PL only implemented for gnls opt.method")
    
    k=dil.x/dil.y
    stopifnot(length(xvar)==length(yvar))
    n=length(xvar)    
    xvar=as.vector(xvar); yvar=as.vector(yvar) # strip off potential attributes; otherwise it may cause trouble with gnls
    dat=data.frame(readout.x=unname(xvar), readout.y=unname(yvar))    
    
    res=list(dilution.ratio=k, dilution.x=dil.x, dilution.y=dil.y, xvar=xvar, yvar=yvar)
    class(res)=c("prc",class(res))
    
    #### Step 1: find init by fitting the naive model assuming no measurement error in variable
    
    if (is.null(init)) {
        init=c("c"=exp(min(xvar,yvar))*0.8, "d"=exp(max(xvar,yvar)), "b"=-1) # *0.8 is to make c smaller than the smallest value of data, assuming readouts can only be positive
        if (!f.is.1) init=c(init,f=1)
        if (verbose) {cat("init:", init, "\n")}    
        if(init.method=="gnls") {
            # note that one of the good things about gnls is that NaN is removed, so c and d estimate can be more reasonable
            if (!f.is.1) {
                formula.gnls = as.formula(  "readout.y ~ log(c+(d-c)/(1+"%+%k%+%"^b*(((d-c)/(exp(readout.x)-c))^(1/f)-1))^f)"  ) 
            } else formula.gnls = as.formula(  "readout.y ~ log(c+(d-c)/(1+"%+%k%+%"^b*(((d-c)/(exp(readout.x)-c))-1)))"  ) 
            # suppress warnings from the following gnls
            fit.1=try(suppressWarnings(gnls(formula.gnls, data=dat, start=init, control=gnlsControl(
                nlsTol=1e-1,  # nlsTol seems to be important, if set to 0.01, then often does not converge
                tolerance=1e-4, 
                # msTol=1e-1, minScale=1e-1, .relStep=1e-7,
                returnObject=TRUE, # allow the return of the fit when the max iter is reached
                maxIter=5000, nlsMaxIter=50, opt="nlminb", msVerbose=T))), silent=FALSE)   
            if(inherits(fit.1, "try-error")) {
                res$error="cannot use naive method (gnls) to find initial valu"
                return (res)
            }
            theta=coef(fit.1)
            if (f.is.1) theta=c(theta,f=1)            
        } else if (init.method=="optim") {
            optim.out = optim(par=init, 
                  fn = function(theta,...) sum(m.0(theta[1],theta[2],theta[3],theta[4],...)), 
                  gr = function(theta,...) colSums(attr(m.0(theta[1],theta[2],theta[3],theta[4],...), "gradient")), 
                  (dat$readout.x), (dat$readout.y), k, 
                  method="BFGS", control = list(trace=0), hessian = F)
            theta=optim.out$par            
        }            
    } else {
        theta=init
    }
    
    if (theta["c"]<0) {
        res$error="the initial estimate of theta has negative c"
        cat(res$error, "\n", file=stderr())
        return (res)        
    }
    
    if (method=="naive") {
        res$coefficients=theta
        res$sigma.sq.naive= mean(m.0(theta["c"], theta["d"], theta["b"], theta["f"], xvar, yvar, k))
        # find rvar and compute orthogonal distance to the curve
        c=theta["c"]; d=theta["d"]; b=theta["b"]; f=theta["f"]
        rvar = numeric(n)
        for (i in 1:n) {
            optimize.out = optimize(m.deriv.r, interval = c(log(c), log(d)), c,d,b,f, xvar[i], yvar[i], k)                           
            rvar[i]=optimize.out$minimum
        }
        res$rvar=rvar
        tmp.3 = m.f (c,d,b,f, rvar, xvar, yvar, k)
        mvar = as.vector(tmp.3)
        res$sigma.sq = mean(mvar) 
       
        return (res)
    }
    
    
    #### Step 2 & 3: loop
    
    if (verbose) cat("LS fit:", theta, "\n")
    rvar = numeric(n)
    # to use gnls, we need to reparameterize to use s.halfk.rvar, which is prc evaluated at rvar with k=sqrt(k)
    s.sqrt.k.rvar = numeric(n)
    iterations=0
    
    
    while(TRUE) {
    
        iterations=iterations+1
        
        dis=numeric(n)
        # Step 2: estimate invidivual mean at a dilution between the two dilutions
        for (i in 1:n) {
            # use optimize
            # interval = c(log(theta["c"]), log(theta["d"])) # this leads to bad convergence behavior
            interval = c(log(theta["c"])-10, log(theta["d"])+10) # a very relaxed constraint while looping
            optimize.out = suppressWarnings(optimize(m.deriv.r, interval = interval, theta["c"], theta["d"], theta["b"], theta["f"], xvar[i], yvar[i], k))
            rvar[i]=optimize.out$minimum
            dis[i]=optimize.out$objective
            
#            optim.out = optim(
#                  par=log(sqrt(theta["c"] * theta["d"])), 
#                  fn = function(x,...) (m.deriv.r(x,...)), 
#                  gr = function(x,...) (attr(m.deriv.r(x,...), "gradient")), 
#                  theta["c"], theta["d"], theta["b"], theta["f"], xvar[i], yvar[i], k, 
#                  method="BFGS", control = list(trace=0), hessian = F)
#            rvar[i]=optim.out$par
        }
        s.sqrt.k.rvar = four_pl_prc (theta["c"], theta["d"], theta["b"], theta["f"], rvar, k^.5) # an easy to make mistake is to use k*.5 instead of k^.5
        
        # Step 3: estimate theta
        if(opt.method=="gnls") {
            dat.stacked=data.frame(readout=c(unname(xvar),unname(yvar)), x=rep(s.sqrt.k.rvar,2), k=rep(c(k^(-1/2),k^(1/2)),each=n))
            if (!f.is.1) {
                formula.gnls = as.formula(  "(readout) ~ log(c+(d-c)/(1+k^b*(((d-c)/(exp(x)-c))^(1/f)-1))^f)"  ) 
                start=theta
            } else {
                formula.gnls = as.formula(  "(readout) ~ log(c+(d-c)/(1+k^b*(((d-c)/(exp(x)-c))-1)))"  ) # 3P
                start=theta[names(theta)!="f"]
            }
            # use suppressWarnings to not print warning msgs like "step halving factor reduced below minimum in NLS step"
            fit.1=try(suppressWarnings(gnls(formula.gnls, data=dat.stacked, start=start, control=gnlsControl(
                nlsTol=1e-1,  # nlsTol seems to be important, if set to 0.01, then often does not converge
                tolerance=1e-4, 
                # msTol=1e-1, minScale=1e-1, .relStep=1e-7,
                returnObject=TRUE, # allow the return of the fit when the max iter is reached
                maxIter=5000, nlsMaxIter=50, opt="nlminb", msVerbose=T))), silent=TRUE)       
            if (inherits(fit.1, "try-error")) {
                res$error="Failure in step 3 gnls call"
                cat(res$error, "\n", file=stderr())
                break;
            }
            new.theta=coef(fit.1)
            if (f.is.1) new.theta=c(new.theta,f=1)
            
        } else if (opt.method=="optim") {
            optim.out = optim(
                  par=theta, 
                  fn = function(theta,...) sum(m.deriv.theta(theta[1],theta[2],theta[3],theta[4],...)), 
                  gr = function(theta,...) colSums(attr(m.deriv.theta(theta[1],theta[2],theta[3],theta[4],...), "gradient")), 
                  (rvar), xvar, yvar, k, 
                  method="BFGS", control = list(trace=0), hessian = F)
            new.theta = optim.out$par 
        
        }
            
        if (verbose) {
            cat("Iter "%+%iterations%+%":", new.theta)
            if (opt.method=="gnls") cat(" logLik:", fit.1$logLik)
            cat("\n")            
        }
        if (max(abs(1 - new.theta/theta)) < reltol) {
            if (verbose) cat("converged\n")
            theta = new.theta # update theta 
            break;
        } else if (iterations>max.iter) {
            if (verbose) cat("max iter reached\n")
            theta = new.theta # update theta 
            break;
        } else {
            theta = new.theta # update theta 
        }        
            
    } # end while loop
    if (!is.null(res$error)) {
        return (res)
    }
    
    c=theta["c"]; d=theta["d"]; b=theta["b"]; f=theta["f"]
    
    # enforce the constraint that rvar are between log(c) and log(d)
    if (verbose) cat("number of rvar needed to be adjusted:", sum(rvar<log(c) & rvar>log(d)), "\n")
    for (i in 1:n) {
        if (rvar[i]<log(c) | rvar[i]>log(d)) { # changed from & to | on Oct 25
            optimize.out = optimize(m.deriv.r, interval = c(log(c), log(d)), c,d,b,f, xvar[i], yvar[i], k)                           
            rvar[i]=optimize.out$minimum
        }
    }
    
    tmp.1 = s.dot.f (c,d,b,f, rvar, k)
    tmp.2 = s.f (c,d,b,f, rvar, k)
    svar = as.vector(tmp.2)
    tmp.3 = m.f (c,d,b,f, rvar, xvar, yvar, k)
    mvar = as.vector(tmp.3)
    
    res$xhat=rvar
    res$yhat=svar
    res$rvar=rvar # so that both naive method and this method have rvar which can be used in predict.prc, in naive method, there is no xhat, just rvar
    res$coefficients=theta # this is here just to have certain order in the fields of res
    res$sigma.sq = mean(mvar) # 2 is asymptotic bias factor    
    
    #### compute asymptotic variance of theta_hat and xhat and yhat
    
    if (model=="4P") {  
        sdot = as.vector(tmp.1)
        sddt = drop(attr(tmp.1, "gradient")[,"r"])
        sddd = drop(attr(tmp.1, "hessian")[,"r","r"])
        dsdot_r_theta.dtheta = attr(tmp.1,"gradient")[,cdbf]
        dsdot_r_theta.dtheta.dtheta = attr(tmp.1,"hessian")[,cdbf,cdbf]
        dsddt_r_theta.dtheta = attr(tmp.1,"hessian")[,"r",cdbf]
        
        ds_r_theta.dtheta = attr(tmp.2,"gradient")[,cdbf]
        ds_r_theta.dtheta.dtheta = attr(tmp.2,"hessian")[,cdbf,cdbf]
        
        A = 1+sdot^2-(yvar-svar)*sddt
        dr.dtheta = A^(-1) * ( (yvar-svar)*dsdot_r_theta.dtheta - sdot*ds_r_theta.dtheta )
            
        ds.dtheta = ds_r_theta.dtheta + sdot * dr.dtheta
        dsdot.dtheta = dsdot_r_theta.dtheta + sddt * dr.dtheta
        dsddt.dtheta = dsddt_r_theta.dtheta + sddd * dr.dtheta
        
        d.ds.dtheta.dtheta = ds_r_theta.dtheta.dtheta + dsdot_r_theta.dtheta %.% dr.dtheta
        d.dsdot.dtheta.dtheta = dsdot_r_theta.dtheta.dtheta + dsddt_r_theta.dtheta %.% dr.dtheta
        
        dr.dtheta.dtheta = -A^(-2) * ( ((yvar-svar)*dsdot_r_theta.dtheta - sdot*ds_r_theta.dtheta) %.% (2*sdot*dsdot.dtheta + ds.dtheta*sddt - (yvar-svar)*dsddt.dtheta) ) + 
            A^(-1) * ( -dsdot_r_theta.dtheta %.% ds.dtheta + (yvar-svar)*d.dsdot.dtheta.dtheta - ds_r_theta.dtheta %.% dsdot.dtheta - sdot * d.ds.dtheta.dtheta )
        
        dm_r_theta.dtheta = attr(tmp.3,"gradient")[,cdbf]
        dm_r_theta.dr = attr(tmp.3,"gradient")[,"r"]    
        dm_r_theta.dtheta.dtheta = attr(tmp.3,"hessian")[,cdbf,cdbf]
        dm_r_theta.dr.dtheta = attr(tmp.3,"hessian")[,"r",cdbf]
        dm_r_theta.dr.dr = attr(tmp.3,"hessian")[,"r","r"]
        
        dm.dtheta = dm_r_theta.dtheta + dm_r_theta.dr * dr.dtheta
        
        dm.dtheta.dtheta = dm_r_theta.dtheta.dtheta + dm_r_theta.dr.dtheta %.% dr.dtheta + dr.dtheta %.% dm_r_theta.dr.dtheta +
            dm_r_theta.dr.dr * (dr.dtheta %.% dr.dtheta) + dm_r_theta.dr * dr.dtheta.dtheta
        
        res$vxhat = res$sigma.sq/(1+sdot^2)
        res$vyhat = res$sigma.sq/(1+sdot^(-2))
            
        V.inv = try(solve(colMeans(dm.dtheta.dtheta)), silent=TRUE)
        if (inherits(V.inv, "try-error")) {
            if (verbose) cat("Fails to estimate covariance matrix\n")
            res$Sigma.hat=NULL
        } else {
            Sigma.hat = 1/n * (V.inv %*% colMeans (dm.dtheta %.% dm.dtheta) %*% V.inv) # 1/n is to get variance of theta^hat
            res$Sigma.hat = Sigma.hat
        }
        
#        # don't delete may be useful later
#        # compute the variance of r_hat, first try, not a good estimator
#        B = (1+sdot^2-{yvar-svar}*sddt)^(-1)
#        dr.dx = B * {1+yvar*sdot}
#        dr.dy = B * sdot
#        res$var.r = (dr.dx^2 + dr.dy^2) * res$sigma.sq #
#        res$var.r = special.mat.f.1 (dr.dtheta, Sigma.hat, dr.dtheta)        
        
        # compute the variance of r_hat, second try
                    
    
    } # end compute asymptotic variance
    
    
    # return object
    res
}    

lines.prc=function(x, col=1, ...) plot(x, type="l", add=TRUE, diag.line=FALSE, lcol=col, ...)
plot.prc=function(x, type=c("b","l","p"), add=FALSE, diag.line=TRUE, lcol=2, pcol=1, log.axis=TRUE, xlab=NULL, ylab=NULL, lwd=2, xlim=NULL, ylim=NULL, ...) {
    
    type <- match.arg(type)    
    coef.=coef(x)
    c=coef.["c"]; d=coef.["d"]; f=coef.["f"]; b=coef.["b"]
    xx=exp(seq(log(c),log(d),length=1e3))
    
    if(log.axis) {
        transf = log; log.axis=""
        if (is.null(xlab)) xlab="x"
        if (is.null(ylab)) ylab="y"
    } else {
        if (is.null(xlab)) xlab="x"
        if (is.null(ylab)) ylab="y"
        transf = I; log.axis="xy"
    }
    
    if (!add) {
        if (is.null(xlim)) {
            if (!is.null(x$var)) {
                xlim=transf(range(exp(x$xvar),exp(x$yvar),c,d)); ylim=xlim
            } else {
                xlim=transf(c(c,d)); ylim=xlim
            }
            
        }
        if (type=="b") {
            plot(transf(exp(x$xvar)),transf(exp(x$yvar)), col=pcol, xlim=xlim, ylim=ylim, log=log.axis, xlab=xlab, ylab=ylab, ...)
            lines(transf(xx), transf(exp(four_pl_prc(c, d, b, f, log(xx), x$dilution.ratio))), col=lcol, lwd=lwd, ...)
        } else if (type=="p") {
            plot(transf(exp(x$xvar)),transf(exp(x$yvar)), col=pcol, xlim=xlim, ylim=ylim, log=log.axis, xlab=xlab, ylab=ylab, ...)
        } else if (type=="l") {
            plot(transf(xx), transf(exp(four_pl_prc(c, d, b, f, log(xx), x$dilution.ratio))), type="l", col=lcol, xlim=xlim, ylim=ylim, log=log.axis, xlab=xlab, ylab=ylab, lwd=lwd, ...)
        }        
    
    } else {
        if (type=="b") {
            points(transf(exp(x$xvar)), transf(exp(x$yvar)), col=pcol, ...)
            lines(transf(xx), transf(exp(four_pl_prc(c, d, b, f, log(xx), x$dilution.ratio))), col=lcol, lwd=lwd, ...)
        } else if (type=="p") {
            points(transf(exp(x$xvar)), transf(exp(x$yvar)), col=pcol, ...)
        } else if (type=="l") {
            lines(transf(xx), transf(exp(four_pl_prc(c, d, b, f, log(xx), x$dilution.ratio))), col=lcol, lwd=lwd, ...)
        }
    
    }
    
    if (diag.line) abline(0,1)
    
}

coef.prc=function(object, ...) {
    ret=object$coefficients
    c(
        c=unname(ifelse(is.na(ret["c"]), exp(ret["logc"]), ret["c"])),
        d=unname(ifelse(is.na(ret["d"]), exp(ret["logd"]), ret["d"])),
        ret["b"],
        ret["f"]
    )
    
}

print.prc=function(x, ...) {
    for (a in names(x)) {
        cat(a,"\n")
        
        if (a %in% c("xvar", "yvar", "xhat", "yhat", "vxhat", "vyhat", "rvar", "p", "support")) {
            print("vector of length "%+%length(x[[a]])%+%" ...", quote=FALSE)
        
        } else if (a == "A") {
            print("matrix of dim n by K", quote=FALSE)
        
        } else {
            print(x[[a]], quote=TRUE)
        
        }
        cat("\n")
    }
    
}

predict.prc=function(object, new.dilution, xvar=NULL, dil.x=NULL, ret.sd=FALSE, ...) {
    if (is.null(xvar) & is.null(dil.x)) {
#        if (max(c(new.dilution,object$dilution.x))/min(c(new.dilution,object$dilution.x)) > max(c(new.dilution,object$dilution.y))/min(c(new.dilution,object$dilution.y))) {
#            old.dilution=object$dilution.y
#            oldvar=object$yvar
#        } else {
#            old.dilution=object$dilution.x
#            oldvar=object$xvar
#        }        
        old.dilution=object$dilution.x
        oldvar=object$rvar
    } else if (!is.null(xvar) & !is.null(dil.x)) {
        old.dilution=dil.x
        oldvar=xvar
    } else stop("the input should have both xvar and dil.x or neither")
    
    k=old.dilution/new.dilution
    coef.=coef(object)
    c=coef.["c"]; d=coef.["d"]; f=coef.["f"]; b=coef.["b"]
    res=four_pl_prc(c, d, b, f, oldvar, k)
    res
}



# expressions and deriv3 functions
# x, y and r are all on the log scale
    
m.0.expr <- expression( (y - log(c+(d-c)/(1+(((d-c)/(exp(x)-c))^(1/f)-1)*k^b)^f)) ^2 )
m.0=deriv3(m.0.expr, c("c","d","b","f"), c("c","d","b","f", "x","y","k"))
    
m.expr <- expression( (y - log(c+(d-c)/(1+(((d-c)/(exp(r)-c))^(1/f)-1)*k^b)^f)) ^2  +  (x - r) ^2 )
#m.f=function(c,d,b,f,r,x,y,k)  (y - log(c+(d-c)/(1+(((d-c)/(exp(r)-c))^(1/f)-1)*k^b)^f))^2  +  (x - r)^2 
m.f=deriv3(m.expr, c("c","d","b","f", "r"), c("c","d","b","f", "r","x","y","k")) # gradient is used in computing variance
m.deriv.r=deriv3(m.expr, c("r"), c("r", "c","d","b","f", "x","y","k"))
m.deriv.theta=deriv3(m.expr, c("c","d","b","f"), c("c","d","b","f", "r","x","y","k"))
    
s.expr <- expression( log(c+(d-c)/(1+(((d-c)/(exp(r)-c))^(1/f)-1)*k^b)^f) )
s.f=deriv3(s.expr, c("c","d","b","f", "r"), c("c","d","b","f", "r", "k")) #  not exported. Within the pkg, it is nice to have a function with as simple a name as s.

#four_pl_prc = function(c,d,b,f, xx, k) log(c+(d-c)/(1+(((d-c)/(exp(xx)-c))^(1/f)-1)*k^b)^f) 
# the .Call version is twice as fast as the R version
four_pl_prc = function(c,d,b,f, xx, k) {
    .Call("compute_four_pl_prc", c,d,b,f, xx, k)
    #.Call("compute_four_pl_prc", .as.double(c),.as.double(d),.as.double(b),.as.double(f), .as.double(xx), k) # very slow
}

s.dot.expr <- expression( k^b * exp(r) * ( 1/(exp(r)-c) * (d-c)/(1+(((d-c)/(exp(r)-c))^(1/f)-1)*k^b)^f )^(1/f+1) / (c+ (d-c)/(1+(((d-c)/(exp(r)-c))^(1/f)-1)*k^b)^f) )
# exported function
s.dot.f=deriv3(s.dot.expr, c("c","d","b","f", "r"), c("c","d","b","f", "r", "k")) 
    
cdbf=c("c","d","b","f")

# these can use C implementation
# A is a matrix n x a, B is a matrix n x b, ret an array n x a x b
"%.%" <- function (A, B) {
    n=nrow(A)
    out=array(dim=c(n, ncol(A), ncol(B)))
    for (i in 1:n) {
        out[i,,]=A[i,] %o% B[i,]
    }
    out
}
# A is a matrix n x a, B is a matrix a x a, C is a matrix n x a, return a vector of length n 
special.mat.f.1 <- function (A, B, C) {
    n=nrow(A)
    out=numeric (n)
    for (i in 1:n) {
        out[i]=A[i,,drop=FALSE] %*% B %*% t(C[i,,drop=FALSE])
    }
    out
}
