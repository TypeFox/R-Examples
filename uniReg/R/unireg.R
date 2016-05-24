negloglikFREQ <- function(lambda,prior,tBSB,tySB,sigmaest,tbetaV,tbetaVbeta,Dtilde,beta0,Vinv,rangV,At,constr){
################################################################################
## Calculates the value of the negative (restricted) log-likelihood of lambda. #
################################################################################
  if(prior(exp(lambda))==0){negloglik <- .Machine$double.xmax
        }else{               
            Einv <- tBSB + exp(lambda)*Vinv
            E <- solve(Einv)
            e <- as.numeric((tySB + exp(lambda)*tbetaV)%*%E) 
            if(constr!="none"){
                Vinvtilde <- Dtilde + exp(lambda)*Vinv
                prob2 <- pmvnorm(lower=0, mean=as.numeric(At%*%e), sigma=At%*%E%*%t(At))[[1]]
                term2 <- -log(prob2)
                prob3 <- pmvnorm(lower=0, mean=as.numeric(At%*%beta0), sigma=At%*%solve(Vinvtilde)%*%t(At))[[1]]
                term3 <- log(prob3)
            }else{term2 <- 0
                term3 <- 0
            }
            term1 <- 0.5*log(det(mean(sigmaest)*Einv)) -log(prior(exp(lambda))) - (rangV/2)*lambda - 0.5*t(e)%*%Einv%*%e + 0.5*exp(lambda)*tbetaVbeta #
            
            negloglik <- term1+term2+term3
            
            if(term3==-Inf){negloglik <- .Machine$double.xmax}
            if(negloglik==Inf){negloglik <- .Machine$double.xmax}
            if(negloglik==-Inf){negloglik <- -.Machine$double.xmax}}
    return(negloglik)
}


unisplinem <- function(x,y,w,tBSB,tySB,sigmaest,tbetaV,tbetaVbeta,Dtilde,rangV,B,p,prior,beta0,Vinv,constr,inverse,penalty,tuning,m){
######################################################################################
# Optimizes the tuning factor lambda and calculates the constrained (with mode m) or #
# unconstrained estimate of the coefficient vector.                                  #
######################################################################################
    At <- inverse*unimat(p,m)          # At is the (transposed) matrix of constraints on the beta coefficients if the mode is m
    
    if(penalty=="none"){lambdaopt <- 0
    }else{
        # if tuning=TRUE, lambda ist optimized with constr=constr, otherwise with constr="none"
        if(!tuning){
            opt <- optimize(negloglikFREQ,interval=c(3,10),prior=prior,tBSB=tBSB,tySB=tySB,sigmaest=sigmaest,tbetaV=tbetaV,tbetaVbeta=tbetaVbeta,Dtilde=Dtilde,beta0=beta0,Vinv=Vinv,rangV=rangV,At=At,constr="none")
            lambdaopt <- exp(opt$minimum)       
        }else{
            opt <- optimize(negloglikFREQ,interval=c(3,10),prior=prior,tBSB=tBSB,tySB=tySB,sigmaest=sigmaest,tbetaV=tbetaV,tbetaVbeta=tbetaVbeta,Dtilde=Dtilde,beta0=beta0,Vinv=Vinv,rangV=rangV,At=At,constr=constr)
            lambdaopt <- exp(opt$minimum)
        }
    }
        
    # constructing E_{lambda}^{-1}, E_{lambda} and e_{lambda} with the optimal lambda    
    Einv <- tBSB+ lambdaopt*Vinv
    tySBlBetaV <- as.numeric(tySB + lambdaopt*tbetaV)
    e <- solve(Einv,tySBlBetaV,system=Einv)
    # constrained coefficients are estimated with solve.QP and constraint matrix At; otherwise the vector e is the unconstrained solution 
    if(constr!="none"){coeff <-solve.QP(Dmat=Einv, dvec=tySBlBetaV, Amat=t(At))$solution
        }else{coeff <- e}
    
    fit <- as.numeric(B%*%coeff)
    if(!is.null(w)){value <- sum(w*(y - fit)^2)}else{value <- sum((y - fit)^2)}
    
    return(list(coef=coeff,fitted.values=fit,lambdaopt=lambdaopt,value=value))
}

unireg <- function(x, y, w=NULL, sigma=NULL, a=min(x), b=max(x), g=10, k=3, constr=c("unimodal","none","invuni","isotonic","antitonic"), 
    penalty=c("diff", "none", "sigEmax", "self", "diag"), Vinv=NULL, beta0=NULL, coinc=NULL, tuning=TRUE, abstol=0.01,vari=5,ordpen=2){
#######################################################
# Applies the specified spline regression to x and y. #
#######################################################
    
    if(!is.vector(x,mode="numeric")){stop("x should be a numeric vector.")}
    if(!is.vector(y,mode="numeric")){stop("y should be a numeric vector.")}
    if(!is.null(w) && !is.vector(w,mode="numeric")){stop("w should be NULL or a numeric vector.")}
    if(!is.null(sigma) && !is.vector(sigma,mode="numeric") && !is.matrix(sigma)){stop("sigma should be NULL or a numeric vector or matrix.")}
    if(!(is.vector(a,mode="numeric") && length(a)==1 && is.finite(a))){stop("a should be a finite numeric vector of length 1.")}
    if(!(is.vector(b,mode="numeric") && length(b)==1 && is.finite(b))){stop("b should be a finite numeric vector of length 1.")}
    if(!(g%%1==0 && g>=0 && is.finite(g))){stop("g should be a finite whole number >=0.")}
    if(!(k%%1==0 && k>=0 && is.finite(k))){stop("k should be a finite whole number >=0.")}
    constr <- match.arg(constr)
    penalty <- match.arg(penalty)
    if(!is.null(Vinv) && !(is.matrix(Vinv) && isTRUE(all.equal(dim(Vinv),c(g+k+1,g+k+1))))){stop("Vinv should be NULL or a (g+k+1)x(g+k+1) matrix.")}
    if(!is.null(beta0) && !is.vector(beta0,mode="numeric")){stop("beta0 should be NULL or a numeric vector.")}
    if(!is.null(coinc) && !is.logical(coinc)){stop("coinc should be NULL, TRUE or FALSE.")}
    if(!is.logical(tuning)){stop("tuning should be TRUE or FALSE.")}
    if(!(is.vector(abstol,mode="numeric") && length(abstol)==1 && is.finite(abstol)) && !is.null(abstol)){stop("abstol should be NULL or a finite numeric vector of length 1.")}
    
    if(!is.null(sigma) && !(is.vector(sigma,mode="numeric") && (length(sigma)==1||sum(sigma==sigma[1])==length(sigma))) && !is.null(abstol)){stop("Cannot estimate sigma in case of heteroscedasticicy. Set abstol to NULL.")}
    if(is.null(sigma) && is.null(abstol)){stop("Either sigma or abstol have to be non-NULL.")}
    if(is.null(sigma) && !is.null(abstol) && !all(table(x)>=2)){stop("Provide a startvalue for sigma or at least 2 repeated measurements for each x-value.")}
    if(b<=a){stop("[a,b] is not a proper interval.")}
    if(penalty=="none" && (g+k+1)>length(unique(x))){warning("Parameters not estimable. Reduce g+k or increase number ob observation points.")}
    
    if(penalty=="self" && (is.null(Vinv) || is.null(beta0))){warning("Vinv and beta0 have to be specified, if penalty='self'.")}
    if(penalty!="self" && (!is.null(Vinv) || !is.null(beta0))){warning("Vinv and beta0 will be ignored, if penalty!='self'.")}
    if(penalty!="self" && is.logical(coinc)){warning("coinc has no influence, if penalty!='self'.")}
    if(penalty=="self" && is.null(coinc)){warning("coinc has to be TRUE or FALSE, if penalty='self'.")}

    if(constr=="none" && !tuning){warning("tuning has no influence, if constr='none'.")}
    if(penalty=="none" && !tuning){warning("tuning has no influence, if penalty='none'.")}

    y <- y[order(x)]
    x <- sort(x)
    n <- length(x)    
    dose <- unique(x)
    nod <- length(dose)
    nd <- as.numeric(table(x)) # number of observations at each dose
    p <- g+k+1
 
    inverse <- 1
    if(constr=="invuni"){
        constr <- "unimodal"
        inverse <- -1
    }
    if(constr=="antitonic"){
        constr <- "isotonic"
        inverse <- -1
    }
    
    if(!is.null(w)){
        if(length(w)!=n){stop("w should be a vector of length n.")
            }else{w <- n*w/sum(w)}
        }else{w <- rep(1,n)}
       
    if(is.null(sigma)){
        wssd <- numeric()
        for(d in seq_along(dose)){
            yd <- y[x==dose[d]]
            wd <- w[x==dose[d]]
            wssd[d] <- sum(wd*(yd-weighted.mean(yd,wd))^2)
        }
        varest <- sum(w)/(sum(w)^2-sum(w^2))*sum(wssd)
    }else{varest <- sigma}

    if(penalty!="self"){
        if(penalty=="none"){coinc <- TRUE; ordpen=2; Dtilde <- diag(0,p)
        }else if(penalty=="diff"){coinc <- FALSE; ordpen=ordpen; Dtilde <- diag(1/vari,p)
        }else if(penalty=="sigEmax"){coinc <- TRUE; ordpen=1; Dtilde <- diag(1/vari,p)
        }else if(penalty=="diag"){coinc <- TRUE; ordpen=0; Dtilde <- diag(0,p)
        }
        # penalty matrix
        Dq <- diag(g+k+1)        
        if(ordpen>=1){Dq <- diff(Dq,difference=ordpen)}
        Vinv <- t(Dq)%*%Dq
    }else{Dtilde <- diag(0,p)}
    
    # determination of knots depending on the interval [a,b], number g of inner knots and degree k of the spline
    # if coinc=T there are k coincident knots at the boundaries a and b
    knotseq <- equiknots(a,b,g,k,coinc)

    # B-spline design matrix
    B <- spline.des(knots=knotseq, x=x, ord = k+1, derivs=rep(0,length(x)), outer.ok = T)$design
    
    yold <- y
    scalfactor <- 0.5*(max(yold) - min(yold))
    shiftfactor <- min(yold) + scalfactor
    y <- (yold-shiftfactor)/scalfactor
    sigmaest <- varest/scalfactor^2
    
    if(penalty!="self"){    
        # calculating beta0
        # if the penalty is a sigEmax one, we fit a sigmoidEmax model and evaluate it at the knot averages
        if(penalty=="sigEmax"){
            bounds <- defBnds(mD = 8)$sigEmax      
            sigE <- fitMod(x, y, model="sigEmax", bnds=bounds)
            knotloc <- knotave(knotseq,d=k)
            beta0 <- predict(sigE, newdata=data.frame(x = knotloc))
        #}else if(penalty=="diag"){
        #    lin <- lm(y ~ x)
        #    knotloc <- knotave(knotseq,d=k)
        #    beta0 <- predict(lin, newdata=data.frame(x = knotloc))
        }else{beta0 <- rep(0,dim(B)[2])} # in all other cases the penalty vector beta0 contains only zeros
    }else{beta0 <- (beta0-shiftfactor)/scalfactor}
    
    prior <- function(x){return(1)}
      
    tbetaVbeta <- t(beta0)%*%Vinv%*%beta0
    tbetaV <- t(beta0)%*%Vinv
    rangV <- qr(Vinv)$rank 
    
    variter <- 0
    repeat{
        tBSB <- t(B)%*%diag(w/sigmaest)%*%B
        tySB <- t(y)%*%diag(w/sigmaest)%*%B
        
        if(constr=="unimodal"){    
            vals <- numeric(p)
            regs <- list()
            for (i in 1:p){
                regs[[i]] <- unisplinem(x,y,w,tBSB=tBSB,tySB=tySB,sigmaest=sigmaest,tbetaV=tbetaV,tbetaVbeta=tbetaVbeta,Dtilde=Dtilde,rangV=rangV,B=B,p=p,prior=prior,beta0=beta0,Vinv=Vinv,constr=constr,inverse=inverse,penalty=penalty,tuning=tuning,i)
                vals[i] <- sum(w*(y - regs[[i]]$fitted.values)^2)
            }
            mini <- which.min(vals)    
            coef.unimod <- regs[[mini]]$coef
            fit.unimod <- regs[[mini]]$fitted.values
            lambdaopt <- regs[[mini]]$lambdaopt   
        }else{ # constr = isotonic or none
            erg <- unisplinem(x,y,w,tBSB=tBSB,tySB=tySB,sigmaest=sigmaest,tbetaV=tbetaV,tbetaVbeta=tbetaVbeta,Dtilde=Dtilde,rangV=rangV,B=B,p=p,prior=prior,beta0=beta0,Vinv=Vinv,constr=constr,inverse=inverse,penalty=penalty,tuning=tuning,p)
            coef.unimod <- erg$coef
            fit.unimod <- erg$fitted.values
            lambdaopt <- erg$lambdaopt
        }
        
        if(!is.null(sigma) && is.null(abstol)){break}
        
        resd <- y-fit.unimod
        varest <- sum(w)/(sum(w)^2-sum(w^2))*sum(w*(resd-weighted.mean(resd,w))^2)
        variter <- variter+1
        if(abs(varest - sigmaest)<abstol || variter>=10){
            sigmaest <- varest
            break}
        sigmaest <- varest
    }
            
    unimod.func <- function(z){
        Bneu <- spline.des(knots=knotseq, z, ord = k+1, derivs=rep(0,length(z)), outer.ok = T)$design
        return(Bneu%*%(scalfactor*coef.unimod+shiftfactor))
    }
    
    return(list(coef=scalfactor*coef.unimod+shiftfactor,x=x,fitted.values=scalfactor*fit.unimod+shiftfactor,unimod.func=unimod.func,lambdaopt=lambdaopt,sigma=sigmaest*scalfactor^2,degree=k,knotsequence=knotseq,g=g,a=a,b=b,variter=variter))
}
