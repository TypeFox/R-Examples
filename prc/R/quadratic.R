# minimize a quadratic mean function

m.quad.0.expr <- expression( (y - (a*x^2+b*x+c)) ^2 )
m.quad.0=deriv3(m.quad.0.expr, c("a","b","c"), c("a","b","c", "x","y"))

m.quad.expr <- expression( (y - (a*r^2+b*r+c)) ^2  +  (x - r) ^2 )
m.quad.f=deriv3(m.quad.expr, c("a","b","c", "r"), c("a","b","c", "r","x","y"))
m.quad.deriv.r=deriv3(m.quad.expr, c("r"), c("r", "a","b","c", "x","y"))
m.quad.deriv.theta=deriv3(m.quad.expr, c("a","b","c"), c("a","b","c", "r","x","y"))


quad.f=function(a,b,c,x) (a*x^2+b*x+c)
    

# init=NULL; reltol=1e-3; opt.method="gnls"; stop.after.init=FALSE; max.iter=50; verbose=TRUE
quadratic.eiv=function(xvar, yvar, init=NULL, reltol=1e-3, opt.method=c("optim"), stop.after.init=FALSE, max.iter=50, verbose=FALSE) {    
    opt.method=match.arg(opt.method)
    if (verbose) {cat("opt.method: ", opt.method, "\n")}
    
    stopifnot(length(xvar)==length(yvar))
    n=length(xvar)    
    xvar=as.vector(xvar); yvar=as.vector(yvar) # strip off potential attributes; otherwise it may cause trouble with gnls
    dat=data.frame(readout.x=unname(xvar), readout.y=unname(yvar))    
    
    res=list(xvar=xvar, yvar=yvar)
    class(res)=c("quad",class(res))
    
    #### Step 1: find init
    
#    if (is.null(init)) init=c("a"=0, "b"=0,"c"=0) # *0.8 is to make c smaller than the smallest value of data, assuming readouts can only be positive
#    if (verbose) {cat("Starting with:", init, "\n")}
#    
#    if(opt.method=="gnls") {
#        # note that one of the good things about gnls is that NaN is removed, so c and d estimate can be more reasonable
#        formula.gnls = as.formula(  "readout.y ~ a*readout.x^2+b*readout.x+c"  ) 
#        fit.1=try(gnls(formula.gnls, data=dat, start=init, control=gnlsControl(
#            nlsTol=1e-1,  # nlsTol seems to be important, if set to 0.01, then often does not converge
#            tolerance=1e-4, 
#            # msTol=1e-1, minScale=1e-1, .relStep=1e-7,
#            maxIter=5000, nlsMaxIter=50, opt="nlminb", msVerbose=T)), silent=FALSE)       
#        theta=coef(fit.1)
#        
#    } else if (opt.method=="optim") {
#        optim.out = optim(par=init, 
#              fn = function(theta,...) sum(m.quad.0(theta[1],theta[2],theta[3],...)), 
#              gr = function(theta,...) colSums(attr(m.quad.0(theta[1],theta[2],theta[3],...), "gradient")), 
#              (dat$readout.x), (dat$readout.y), 
#              method="BFGS", control = list(trace=0), hessian = F)
#        theta=optim.out$par
#        
#    }            
#    if (verbose) cat("Initial theta:", theta, "\n")
#    if (stop.after.init) {
#        if (verbose) cat("stop.after.init\n")
#        res$coefficients=theta
#        res$sigma.sq= mean(m.quad.0(theta["a"], theta["b"], theta["c"], xvar, yvar))
#        return (res)
#    }
    
    if (!is.null(init)) {
        theta=init        
    } else {
        init=c("a"=0, "b"=0,"c"=0) # *0.8 is to make c smaller than the smallest value of data, assuming readouts can only be positive
        if(opt.method=="gnls") {
            # note that one of the good things about gnls is that NaN is removed, so c and d estimate can be more reasonable
            formula.gnls = as.formula(  "readout.y ~ a*readout.x^2+b*readout.x+c"  ) 
            fit.1=try(gnls(formula.gnls, data=dat, start=init, control=gnlsControl(
                nlsTol=1e-1,  # nlsTol seems to be important, if set to 0.01, then often does not converge
                tolerance=1e-4, 
                # msTol=1e-1, minScale=1e-1, .relStep=1e-7,
                maxIter=5000, nlsMaxIter=50, opt="nlminb", msVerbose=T)), silent=FALSE)       
            theta=coef(fit.1)
            
        } else if (opt.method=="optim") {
            optim.out = optim(par=init, 
                  fn = function(theta,...) sum(m.quad.0(theta[1],theta[2],theta[3],...)), 
                  gr = function(theta,...) colSums(attr(m.quad.0(theta[1],theta[2],theta[3],...), "gradient")), 
                  (dat$readout.x), (dat$readout.y), 
                  method="BFGS", control = list(trace=0), hessian = F)
            theta=optim.out$par
            
        }            
        if (verbose) cat("Initial theta:", theta, "\n")
        if (stop.after.init) {
            if (verbose) cat("stop.after.init\n")
            res$coefficients=theta
            res$sigma.sq= mean(m.quad.0(theta["a"], theta["b"], theta["c"], xvar, yvar))
            return (res)
        }
    }
    if (verbose) {cat("Starting with:", init, "\n")}
    
    
    
    
    #### Step 2 & 3: loop
    
    rvar = numeric(n)
    # to use gnls, we need to reparameterize to use s.halfk.rvar, which is quad evaluated at rvar with k=sqrt(k)
    s.sqrt.k.rvar = numeric(n)
    iterations=0
    
    while(TRUE) {
    
        iterations=iterations+1
        
        # Step 2: estimate invidivual mean at a dilution between the two dilutions
        for (i in 1:n) {
            # use optimize
            optimize.out = suppressWarnings(optimize(m.quad.deriv.r, interval=c(-10,10), theta["a"], theta["b"], theta["c"], xvar[i], yvar[i]))
            rvar[i]=optimize.out$minimum            
        }
        
        # Step 3: estimate theta
        if(opt.method=="gnls") {
            
        } else if (opt.method=="optim") {
            optim.out = optim(
                  par=theta, 
                  fn = function(theta,...) sum(m.quad.deriv.theta(theta[1],theta[2],theta[3],...)), 
                  gr = function(theta,...) colSums(attr(m.quad.deriv.theta(theta[1],theta[2],theta[3],...), "gradient")), 
                  rvar, xvar, yvar,
                  method="BFGS", control = list(trace=0), hessian = F)
            new.theta = optim.out$par 
        
        }
            
        if (verbose) {
            cat("Iter "%+%iterations%+%".", "theta:", new.theta)
            if (opt.method=="optim") cat(" val:", optim.out$value)
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
        
    a=theta["a"]; b=theta["b"]; c=theta["c"]    
    tmp.3 = m.quad.f (a,b,c, rvar, xvar, yvar)
    mvar = as.vector(tmp.3)
    
    # return object
    res$coefficients=theta
    res$rvar=rvar
    res$sigma.sq = mean(mvar) # 2 is asymptotic bias factor    
    res
}    

coef.quad=function(object, ...) {
    ret=object$coefficients
    c(ret[1],ret[2],ret[3])
    
}

print.quad=function(x, ...) {
    for (a in names(x)) {
        cat(a,"\n")
        
        if (a %in% c("xvar", "yvar", "rvar")) {
            print("vector of length "%+%length(x[[a]])%+%" ...", quote=FALSE)
        
        } else if (a == "A") {
            print("matrix of dim n by K", quote=FALSE)
        
        } else {
            print(x[[a]], quote=TRUE)
        
        }
        cat("\n")
    }
    
}

lines.quad=function(x, col=1, x.range=NULL, ...) plot(x, type="l", add=TRUE, lcol=col, x.range=x.range, ...)
plot.quad=function(x, type=c("b","l"), add=FALSE, lcol=2, pcol=1, xlab=NULL, ylab=NULL, lwd=2, x.range=NULL, log.axis=TRUE, ...) {
    
    type <- match.arg(type)    
    coef.=coef(x)
    a=coef.["a"]; b=coef.["b"]; c=coef.["c"]
    if (is.null(x.range)) x.range=c(-100,100)
    xx=seq(x.range[1],x.range[2],length=1e3)
    
    transf=I    
    if (!add) {
        if (type=="b") {
            plot(transf((x$xvar)),transf((x$yvar)), col=pcol, xlab=xlab, ylab=ylab, ...)
            lines(transf(xx), transf((quad.f(a,b,c,xx))), col=lcol, lwd=lwd, ...)
        } else if (type=="l") {
            plot(transf(xx), transf((quad.f(a,b,c,xx))), type="l", col=lcol, log=log.axis, xlab=xlab, ylab=ylab, lwd=lwd, ...)
        }        
    
    } else {
        if (type=="b") {
            points(transf((x$xvar)), transf((x$yvar)), col=pcol, ...)
            lines(transf(xx), transf((quad.f(a,b,c,xx))), col=lcol, lwd=lwd, ...)
        } else if (type=="l") {
            lines(transf(xx), transf((quad.f(a,b,c,xx))), col=lcol, lwd=lwd, ...)
        }
    
    }
        
}
