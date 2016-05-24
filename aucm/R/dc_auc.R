smooth.rauc=function(formula, data, ...) auc.dca(formula, data, type="srauc" ,...)
srauc=smooth.rauc

dcsauc=function(formula, data, ...) auc.dca(formula, data, type="dcsauc" ,...)


# dca for auc maximization
auc.dca=function(formula, data, 
    type="srauc",
    kernel="linear", para=NULL,  
    lambda=.1, zeta=.1, b=10, s=1, epsilon=1e-3, # sauc-dca parameters
    method="tron", decomposition=TRUE,
    dca.control =  list(maxit=1e3, abstol=1e-5, coef.init=NULL, lincomb.init=NULL),
    tron.control = list(q=50, maxfev=1e3, gtol=1e-2, frtol=1e-12, K.thresh=1, verbose=0),
    return.K=FALSE, verbose = FALSE
){ 
    
    begin=Sys.time()        
    
    if (!type %in% c("srauc","dcsauc")) stop("type not supported") 
    
    if (type=="srauc" & !decomposition) {
        decomposition=TRUE
        warning("for srauc and tron, non-decomposition is not implemented, changing to decomposition")
    }
    
    if (is.null(dca.control$maxit)) dca.control$maxit=1e3
    if (is.null(dca.control$abstol)) dca.control$abstol=1e-5
    
    if (is.null(tron.control$q)) tron.control$q=50
    if (is.null(tron.control$maxfev)) tron.control$maxfev=1e3
    if (is.null(tron.control$gtol)) tron.control$gtol=1e-2 # setting this to 1e-4 slows things down significantly
    if (is.null(tron.control$frtol)) tron.control$frtol=1e-12
    if (is.null(tron.control$K.thresh)) tron.control$K.thresh=1
    if (is.null(tron.control$verbose)) tron.control$verbose=0
    
    # method can be tron or other methods supported by optim
    method = tolower(substr(method[1],1,1))
    if(!method %in% c('t')){
        method = match(method,c("n","b","c","l","s"),nomatch = NA)
        if(is.na(method))stop("Please specify one of the optimization methods Tron(t/T) Nelder-Mead(n/N), BFGS(b/B), CG(c/C), L-BFGS-B(l/L), SANN=(s/S).")
        method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")[method]
    }
    
    tmp=model.frame(formula, data)
    X1=model.matrix(formula, tmp[tmp[,1]==1,])[,-1,drop=FALSE]
    X2=model.matrix(formula, tmp[tmp[,1]==0,])[,-1,drop=FALSE]
    n.case=nrow(X1)
    n.control=nrow(X2)
    n=n.case+n.control
    n.1.n.2=n.case*n.control
    X=rbind(X1,X2)
    # We stack case on top of control. To access j in the control, we add n.case
    
    p = NULL
    beta.k  = NULL
    theta.k = NULL
    X.diff  = NULL
    K       = NULL
    
    kernel=substr(kernel,1,1)
    if (kernel!="l" && is.null(para)) stop("kernel parameter not set") 
    
    pairs.to.exclude = 0
    if(kernel == 'l'){    
    
        p = ncol(X)
        K = diag(p)
        X.diff=get.X.diff(X1,X2)
        
        if (is.null(dca.control$coef.init)) {
            fit.rlogit=rlogit(formula, data, verbose=verbose)
            if (fit.rlogit$convergence) {
                theta.k=fit.rlogit$coef[-1]
            } else{ 
                theta.k  = double(p) + 1.0
            }                        
        } else {
            theta.k = dca.control$coef.init
        }
        
        if (verbose>=2) {
            cat("Initial values: ", theta.k)
        }
        
    } else {      
    
        p = n
        K=getK(X,kernel,para)
        
        X.diff=K[rep(1:n.case,each=n.control), ,drop=FALSE] - K[rep(n.case+1:n.control,n.case),,drop=FALSE] # it is weird to call it X.diff here, but it stuck
        
        if (!is.null(dca.control$lincomb.init)) {
            theta.k = solve(K, dca.control$lincomb.init)
            myprint(theta.k)
        } else {
            if (is.null(dca.control$coef.init)) {
                theta.k = double(p) + 0.0
            } else {
                theta.k = dca.control$coef.init
            }        
        }        
    
        if (method=="t") {
            # get pairs where K[i,j] is too close to 1, which may bring trouble to tron
            pairs.high = which(K<1 & K >tron.control$K.thresh, arr.ind=TRUE)
            pairs.high = pairs.high[pairs.high[,1]<pairs.high[,2],]
            pairs.high = pairs.high-1 # change to 0-based index for passing to .C
            pairs.to.exclude = nrow(pairs.high)
            pairs.to.exclude.a = as.integer(pairs.high[,1])
            pairs.to.exclude.b = as.integer(pairs.high[,2])
            #if (verbose) print (pairs.high)
        }
    
    }
    
    # cannot pass null to .C
    if (pairs.to.exclude==0) {
        pairs.to.exclude.a=integer(1)
        pairs.to.exclude.b=integer(1)
    }
    
    
    rel.err = Inf
    iter=0
    tron.convergence <- as.integer(0)
    penalized.losses <- losses <- NULL
    
    f.dif.k = X.diff %*% theta.k
    if (type=="srauc") linear.coef.in.dca=ifelse(f.dif.k < -.5*s, -1/s, 0.0) else linear.coef.in.dca=S2.deriv(f.dif.k, zeta, b) # initialize
    
    
    while(TRUE) {
    
        iter=iter+1
        if (iter>dca.control$maxit) break
    
#        if(iter>2) {
#            # plot likelihood surface
#            require(rgl)
#            xx=expand.grid(seq(1,6,length=10), seq(-6,-1,length=10))
#            z=sapply(1:nrow(xx), function (xi) rauc.dca.approx(t(xx[xi,]), X.diff, linear.coef.in.dca, epsilon))
#            print(summary(z))
#            #plot3d(xx[,1], xx[,2], z)
#            print(str(z))
#            par(mfrow=c(1,2))
#            plot(xx[51:60,1], z[51:60], type="l")
#            plot(xx[0:9*10+5,2], z[0:9*10+5], type="l")
#        }
#            
    
        if (method=="t") {
            # use tron and tron-decomposition to find solution
            if (type=="dcsauc") {
                if (decomposition & tron.control$q>p) tron.control$q=p #stop("this cannot be: q>p "%+%tron.control$q%+%" "%+%p)
                param=theta.k * 1 # *1 forces a value-copy instead of reference-copy
                # using if in .C first argument leads to an error in R CMD check
                if(decomposition)
                    optim.res <- .C("dcsauc_tron_decomposition",
                        K, X.diff, linear.coef.in.dca, lambda, as.integer(p), as.integer(n.1.n.2), zeta, b,
                        # tron control variables. maxfev is similar to maxit, gtol is relative tolerance
                        as.integer(tron.control$maxfev), tron.control$gtol, tron.control$frtol, as.integer(tron.control$verbose), as.integer(tron.control$q),
                        pairs.to.exclude, pairs.to.exclude.a, pairs.to.exclude.b, # defined as integers
                        # output variables
                        par=param, value = double(1), convergence=integer(1), 
                        # house keeping variables
                        DUP = TRUE, NAOK = FALSE, PACKAGE = "aucm"
                    )             
                else 
                    optim.res <- .C("dcsauc_tron",
                        K, X.diff, linear.coef.in.dca, lambda, as.integer(p), as.integer(n.1.n.2), zeta, b,
                        # tron control variables. maxfev is similar to maxit, gtol is relative tolerance
                        as.integer(tron.control$maxfev), tron.control$gtol, tron.control$frtol, as.integer(tron.control$verbose), as.integer(tron.control$q),
                        pairs.to.exclude, pairs.to.exclude.a, pairs.to.exclude.b, # defined as integers
                        # output variables
                        par=param, value = double(1), convergence=integer(1), 
                        # house keeping variables
                        DUP = TRUE, NAOK = FALSE, PACKAGE = "aucm"
                    )             
            } else if (type=="srauc") {
                if (decomposition & tron.control$q>p) tron.control$q=p
                param=theta.k * 1 # *1 forces a value-copy instead of reference-copy
                if (decomposition)
                    optim.res <- .C("srauc_tron_decomposition",
                        K, X.diff, linear.coef.in.dca, lambda, as.integer(p), as.integer(n.1.n.2), epsilon,
                        # tron control variables. maxfev is similar to maxit, gtol is relative tolerance
                        as.integer(tron.control$maxfev), tron.control$gtol, tron.control$frtol, as.integer(tron.control$verbose), as.integer(tron.control$q),
                        pairs.to.exclude, pairs.to.exclude.a, pairs.to.exclude.b, # defined as integers
                        # output variables
                        par=param, value = double(1), convergence=integer(1), 
                        # house keeping variables
                        DUP = TRUE, NAOK = FALSE, PACKAGE = "aucm"
                    )             
                else 
                    optim.res <- .C("srauc_tron",
                        K, X.diff, linear.coef.in.dca, lambda, as.integer(p), as.integer(n.1.n.2), epsilon,
                        # tron control variables. maxfev is similar to maxit, gtol is relative tolerance
                        as.integer(tron.control$maxfev), tron.control$gtol, tron.control$frtol, as.integer(tron.control$verbose), as.integer(tron.control$q),
                        pairs.to.exclude, pairs.to.exclude.a, pairs.to.exclude.b, # defined as integers
                        # output variables
                        par=param, value = double(1), convergence=integer(1), 
                        # house keeping variables
                        DUP = TRUE, NAOK = FALSE, PACKAGE = "aucm"
                    )             
            }
        } else {
            # use R optim to find solution
            if (type=="dcsauc") {
                optim.res=optim(par=theta.k, 
                    fn = function(x) sauc.approx.dca(x, X.diff, linear.coef.in.dca, zeta, b)+.5*lambda*c(crossprod(x,K)%*%x), 
                    gr = function(x) sauc.approx.dca.gr(x, X.diff, linear.coef.in.dca, zeta, b) + lambda*c(crossprod(x,K)), 
                    method=method)
            } else if (type=="srauc") {
            
                optim.res=optim(par=theta.k, 
                    fn = function(x) rauc.dca.approx(x, X.diff, linear.coef.in.dca, epsilon) + .5*lambda*c(crossprod(x,K)%*%x), 
                    gr = function(x) rauc.dca.approx.gr(x, X.diff, linear.coef.in.dca, epsilon) + lambda*c(crossprod(x,K)), 
                    method=method)
            }                 
        }
        
        theta.k = optim.res$par * 1
        f.dif.k = X.diff %*% theta.k
        if (type=="dcsauc") linear.coef.in.dca=S2.deriv(f.dif.k, zeta, b) else if (type=="srauc") linear.coef.in.dca=ifelse(f.dif.k < -.5*s, -1/s, 0.0) # recompute
        
        loss = ifelse (type=="dcsauc", mean(dcsauc.f(f.dif.k, zeta, b)), mean(ramp.f(f.dif.k, 1)))
        losses = c(losses, loss)
        penalty   = .5*lambda*c(crossprod(theta.k,K)%*%theta.k)/n.1.n.2
        penalized.losses=c(penalized.losses, loss+penalty)
        
        if(verbose) {
            cat("dca #", iter, 
                ", loss+penalty:",format(round(loss+penalty,6),nsmall=5), 
                ", loss:",loss, 
                ", penalty:",penalty, 
            sep="")            
            if(kernel == 'l') cat("    theta:",concatList(signif(theta.k,5),"|"),"\n") else {
                sm = summary(theta.k); cat("    theta:",paste(names(sm),sm,sep='='),"\n")
            }
            cat("\n");
        }
        
#        # use reltol as stopping criterion
#        if(kernel == 'l'){
#            rel.err = max(abs(optim.res$par - theta.k)/theta.k)
#        } else{
#            rel.err = abs(optim.res$value - objval) / objval
#            objval = optim.res$value
#        }        
#        if (rel.err < dca.control$reltol) break 
    
        # use abstol as stopping criterion
        if (iter>1) {
            delta=penalized.losses[iter]-penalized.losses[iter-1]
            if (abs(delta)<dca.control$abstol) break               
        }
    }
    
    fit=list()
    if (type=="dcsauc") class(fit)=c("dcsauc","auc",class(fit))
    if (type=="srauc") class(fit)=c("srauc","auc",class(fit))
    
    fit$convergence=ifelse (iter>dca.control$maxit, 1, 0)
    fit$converged=ifelse (iter>dca.control$maxit, FALSE, TRUE)
    fit$iterations=iter
    fit$time.elapsed = difftime(Sys.time(),begin,units = 'secs')
    
    fit$coefficients=optim.res$par
    fit$losses = losses 
    fit$penalized.losses = penalized.losses 
    
    fit$formula=formula
    fit$kernel=kernel
    fit$para=para
    fit$X=X
    fit$y=c(rep(1,n.case),rep(0,n.control))
    fit$n.case=n.case
    fit$n.control=n.control
    fit$lambda=lambda
    
    if(kernel == 'l'){
        fit$linear.combination = c(X %*% fit$coefficients)
    } else {
        fit$linear.combination = c(K %*% fit$coefficients)
    }
    fit$train.auc = fast.auc(fit$linear.combination, fit$y)
    
    if (return.K) fit$K=K    
    
    if (type=="dcsauc") {
        tmp=dcsauc.f(f.dif.k, zeta, b)
        #fit$saturation1=try(quantile(tmp, c(.02,.1,.2,.8,.9,.98)))
        fit$saturation=mean(tmp>.99 | tmp<0.01)
    } else {
        fit$saturation=mean(fit$linear.combination>.5*s | fit$linear.combination< -.5*s)
    }
    
    if(verbose){
        cat("saturation: ",fit$saturation,"\n")
        print(fit$time.elapsed)
    }
    fit
    
}



S1=function(x, zeta, b) log(1+exp(-b*x))/zeta
S1.deriv=function(x, zeta, b) -b/(1+exp(b*x))/zeta
S1.dd=function(x, zeta, b) b^2/(1+exp(b*x))/(1+exp(-b*x))/zeta

S2=function(x, zeta, b) log(1+exp(-b*x-zeta))/zeta
S2.deriv=function(x, zeta, b) -b/(1+exp(b*x+zeta))/zeta

S3=function(x, s) log(1+exp(.5-x/s)) # a shifted logistic likelihood function
S3.deriv=function(x, s) -1/s/(1+exp(x/s-0.5))

#S4=function(x, epsilon) 0.5-epsilon*(log(1+exp(0.5/epsilon)) - log(1+exp((0.5-x)/epsilon)))
S4=function(x, epsilon) 0.5-epsilon*(log(1+exp(0.5/epsilon)) - ifelse((0.5-x)/epsilon>700, (0.5-x)/epsilon, log(1+exp((0.5-x)/epsilon)))) # ifelse is to deal with exp(800)=Inf
S4.deriv=function(x, epsilon) -1/(1+exp((x-0.5)/epsilon))
S4.dd = function(x, epsilon) 1/epsilon/(1+exp((x-0.5)/epsilon))/(1+exp((0.5-x)/epsilon))

# a function that approximates sauc
dcsauc.f = function(eta, zeta, b) {
    drop(S1(eta, zeta, b) - S2(eta, zeta, b))
}

# beta and X.diff has to be separate, b/c beta is to be optimized by optim
sauc.approx.dca = function(beta, X.diff, S2.deriv.ij, zeta, b) {
    eta = X.diff %*% beta
    sum(S1(eta, zeta, b) - eta * S2.deriv.ij)
}

sauc.approx.dca.gr = function(beta, X.diff, S2.deriv.ij, zeta, b) {
    eta = X.diff %*% beta
    colSums(drop(S1.deriv(eta, zeta, b) - S2.deriv.ij)*X.diff)
    # t(X.diff) %*% (S1.deriv(eta, zeta, b) - S2.deriv.ij)
}

sauc.approx.dca.hess = function(beta, X.diff, zeta, b) {
    eta = X.diff %*% beta
    tmp = sapply(1:nrow(X.diff), simplify="array", function (i) S1.dd(eta[i], zeta, b) * outer(X.diff[i,], X.diff[i,]) )
    apply(tmp, 1:2, sum)
}

# beta and X.diff has to be separate, b/c beta is to be optimized by optim
rauc.dca.approx = function(beta, X.diff, linear.coef.in.dca, epsilon) {
    eta = X.diff %*% beta    
    sum(S4(eta, epsilon) - eta * linear.coef.in.dca)
}

rauc.dca.approx.gr = function(beta, X.diff, linear.coef.in.dca, epsilon) {
    eta = X.diff %*% beta
    colSums(drop(S4.deriv(eta, epsilon) - linear.coef.in.dca)*X.diff)
}





#    # expressions for use with optim with deriv3
#    eta.s="("%+% concatList("b"%+%1:p%+%"*x"%+%1:p%+%"","+") %+%")"
#    S1.s="log(1+exp(-"%+%eta.s%+%"))/zeta"
#    loss.s=S1.s%+%" + "%+%eta.s%+%"*S2.deriv.ij"
#    loss.f=deriv3(parse(text=loss.s), "b"%+%1:p, c("b"%+%1:p,"x"%+%1:p,"S2.deriv.ij","zeta"))
#    
#        # if we don't want to code the deriviative functions by hand, we can use deriv3, but we need to deal with vector parameters through string manipulation, but this can be a lot slower
#        cri.s="sum(loss.f(" %+% concatList("beta["%+%1:p%+%"]",",") %+% "," %+% concatList("X.diff[,"%+%1:p%+%"]",",") %+% ",S2.deriv.ij,zeta)) + .5*lambda*sum(beta^2)"
#        gr.s="colSums(attr(loss.f(" %+% "" %+% concatList("beta["%+%1:p%+%"]",",") %+% "," %+% concatList("X.diff[,"%+%1:p%+%"]",",") %+% ",S2.deriv.ij,zeta), \"gradient\")) + lambda*beta"                
#        optim.res=optim(par=rep(1,p), fn = function(beta) eval(parse(text=cri.s)), gr = function(beta) eval(parse(text=gr.s)), 
#              method="Nelder-Mead", control = list(), hessian = FALSE)


#       # debugging    
#        beta=beta.k
#        beta=c(-1.7143, 2.8123)    
#        
#        eta = X.diff %*% beta
#        sum(S1(eta, zeta) - eta * S2.deriv.ij - S2(f.dif.k, zeta) + f.dif.k * S2.deriv.ij)
#        sum(S1(eta, zeta) - S2(eta, zeta))
#
#        sum(h1(eta) - eta * h2.deriv(f.dif.k) - h2(f.dif.k) + f.dif.k * h2.deriv(f.dif.k) )
#        sum(h1(eta) - h2(eta))
#
#        S1(eta, zeta) + eta * S2.deriv.ij - S2(f.dif.k, zeta) - f.dif.k * S2.deriv.ij - {S1(eta, zeta) - S2(eta, zeta)} 
#
#        h1(eta) + eta * h2.deriv(f.dif.k) - h2(f.dif.k) - f.dif.k * h2.deriv(f.dif.k) - {h1(eta) - h2(eta)} # diff bt actual and approximation
        
