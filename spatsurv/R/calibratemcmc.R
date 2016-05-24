##' QuadApprox function
##'
##' A function to compute the second derivative of a function (of several real variables) using a quadratic approximation  on a 
##' grid of points defined by the list argRanges. Also returns the local maximum. 
##'
##' @param fun a function
##' @param npts integer number of points in each direction
##' @param argRanges a list of ranges on which to construct the grid for each parameter 
##' @param plot whether to plot the quadratic approximation of the posterior (for two-dimensional parameters only)
##' @param ... other arguments to be passed to fun
##' @return a 2 by 2 matrix containing the curvature at the maximum and the (x,y) value at which the maximum occurs 
##' @export


QuadApprox <- function(fun,npts,argRanges,plot=FALSE,...){

    npar <- length(argRanges)
    vals <- lapply(argRanges,function(x){seq(x[1],x[2],length.out=npts)})
    parn <- paste("x",1:npar,sep="")
    parnames <- parn # later this will be used in the call to lm
    paridx <- as.list(rep(NA,1+2*npar+choose(npar,2))) # this will be used later to find the function maximum via a set of simultaneous equations (intercept, single, squared and mixed terms)
    paridx[(1+1:npar)] <- 1:npar
    gr <- expand.grid(vals)
    gr2 <- gr^2
    parnames <- c(parnames,paste("x",1:npar,".2",sep=""))
    paridx[(npar+2):(2*npar+1)] <- 1:npar
    grcross <- matrix(NA,nrow(gr),choose(npar,2))
    ct <- 1
    for(i in 1:(npar-1)){
        for(j in (i+1):npar){    
            grcross[,ct] <- gr[,i]*gr[,j]
            parnames <- c(parnames,paste(parn[i],parn[j],collapse="",sep=""))
            paridx[[2*npar+1+ct]] <- c(i,j)
            ct <- ct + 1
        }
    }
    partype <- c("intercept",rep("single",npar),rep("squared",npar),rep("mixed",choose(npar,2)))
    dataf <- cbind(gr,gr2,grcross)
    names(dataf) <- parnames
    
    cat("Constructing quadratic approximation to posterior (this can take some time) ...\n")
    dataf$funvals <- apply(gr,1,function(params){fun(params,...)})
    cat("Done.\n")
    
    
    if(plot){
        if(npar==2){
            image.plot(vals[[1]],vals[[2]],matrix(dataf$funvals,npts,npts),main="Function")
        }  
    }
    
    form <- paste("funvals ~",paste(parnames,collapse=" + "))

    mod <- lm(form,data=dataf)
    co <- coefficients(mod)
    
    if(plot){
        if(npar==2){
            image.plot(vals[[1]],vals[[2]],matrix(fitted(mod),npts,npts),main="Quadratic Approximation")
        }  
    }    
    
    # now construct matrix of second derivatives
    sigmainv <- matrix(NA,npar,npar)
    diag(sigmainv) <- 2 * co[which(partype=="squared")] # first the diagonal elements
    idx <- which(partype=="mixed") # now the off diagonals
    ct <- 1
    for(i in 1:(npar-1)){
        for(j in (i+1):npar){    
            sigmainv[i,j] <- co[idx[ct]]
            sigmainv[j,i] <- co[idx[ct]]
            ct <- ct + 1
        }
    }
    
    # lastly, create a system of simultaneous equations, Ax = b, which when solved gives the maximum
    b <- (-1) * matrix(co[which(partype=="single")],npar,1)
    A <- matrix(NA,npar,npar)
    diag(A) <- 2 * co[which(partype=="squared")]
    for(i in 1:(npar-1)){
        for(j in (i+1):npar){
            tst <- sapply(paridx,function(x){any(x==i)&any(x==j)})
            idx <- which(tst)   
            A[i,j] <- co[idx]
            A[j,i] <- co[idx]
        }
    }

    etaest <- as.vector(solve(A)%*%b) # now solve the system of simultaneous equations to get an initial guess for eta 
    
    sigmainv <- fixmatrix(sigmainv)

    return(list(max=etaest,curvature=sigmainv,mod=mod)) 
}




##' fixmatrix function
##'
##' !! THIS FUNCTION IS NOT INTENDED FOR GENERAL USE !!
##'
##' A function to fix up an estimated covariance matrix using a VERY ad-hoc method.
##'
##' @param mat a matrix
##' @return the fixed matrix
##' @export

fixmatrix <- function(mat){
    
    mat <- (-1)*mat # since mat is curvature, the negative *should* have positive eigenvalues
    ev <- eigen(mat)$values
    if(all(ev>0)){
        return((-1)*mat)
    }
    else if(all(ev<0)){
        stop("Estimated covariance matrix for eta has all negative eigenvalues")
    }
    else{
        warning("Something is wrong with the estimated covariance matrix, fixing this using a totally ad-hoc method. This will not affect ergodicity, merely the efficiency of the chain.",immediate.=TRUE)
        cat("Fixing non positive definite covariance matrix for eta ...\n") 
    
        diag(mat) <- abs(diag(mat)) # hmmmm ....        
               
        if(all(dim(mat)==2)){
            fun <- function(x){
                tmp <- mat
                tmp[1,2] <- tmp[2,1] <- mat[1,2] / x
                posev <- abs(ev)
                ev1 <- eigen(tmp)$values
                if(!all(ev1>0)){
                    return(.Machine$double.xmax)
                }
                else{
                    df1 <- (posev[1]-ev1[1])/posev[1]
                    df2 <- (posev[2]-ev1[2])/posev[2]
                    return(df1^2+df2^2)
                }                
            }
            op <- suppressWarnings(try(optimise(fun,interval=c(0,10))))
            if(inherits(op,"try-error")){
                stop("Failed to fix negative definite matix")
            }
            ans <- mat
            ans[1,2] <- ans[2,1] <- mat[1,2] / op$minimum
                       
        }
        else{
            #browser()
            fun1 <- function(pars){
                tmp <- mat
                tmp[lower.tri(tmp)] <- tmp[lower.tri(tmp)] / pars
                tmp[upper.tri(tmp)] <- tmp[upper.tri(tmp)] / pars
                posev <- abs(ev)
                ev1 <- eigen(tmp)$values
                if(!all(ev1>0)){
                    return(.Machine$double.xmax)
                }
                else{
                    dff <- sum(((posev-ev1)/posev)^2)
                    return(dff)
                }                
            }
            op <- suppressWarnings(try(optim(par=rep(1,ncol(mat)),fn=fun1)))
            if(inherits(op,"try-error")){
                stop("Failed to fix negative definite matix")
            }
            
            ans <- mat
            ans[lower.tri(ans)] <- ans[lower.tri(ans)] / op$par
            ans[upper.tri(ans)] <- ans[upper.tri(ans)] / op$par
        }
        

        ct <- nrow(ans)        
        if(!all(eigen(ans)$values>0)){
            while(!all(eigen(ans)$values>0) & ct>0){
                ans[ct,1:(ct-1)] <- 0
                ans[1:(ct-1),ct] <- 0
                ct <- ct - 1
            }    
        }
        
        if(!all(eigen(ans)$values>0)){
            stop("Failed to fix negative definite matix")
        }
        
        ans <- (-1)*ans              
        
        return(ans)    
    } 
}



##' proposalVariance function
##'
##' A function to compute an approximate scaling matrix for the MCMC algorithm. Not intended for general use. 
##'
##' @param X the design matrix, containing covariate information 
##' @param surv an object of class Surv 
##' @param betahat an estimate of beta 
##' @param omegahat an estimate of omega 
##' @param Yhat an estimate of Y 
##' @param priors the priors 
##' @param cov.model the spatial covariance model 
##' @param u a vector of pairwise distances 
##' @param control a list containg various control parameters for the MCMC and post-processing routines 
##' @return an estimate of eta and also an approximate scaling matrix for the MCMC
##' @export

proposalVariance <- function(X,surv,betahat,omegahat,Yhat,priors,cov.model,u,control){
     
    n <- nrow(X)
    lenbeta <- length(betahat)
    lenomega <- length(omegahat)
    leneta <- getleneta(cov.model)
    lenY <- length(Yhat)
    npars <- lenbeta + lenomega + leneta + lenY
    
    sigma <- matrix(0,npars,npars)
    
    # eta
    logpost <- function(eta,surv,X,beta,omega,Y,priors,cov.model,u,control){

        etapars <- cov.model$itrans(eta)
        sigma <- matrix(EvalCov(cov.model=cov.model,u=u,parameters=etapars),n,n)
        cholsigma <- t(chol(sigma))
        cholsigmainv <- solve(cholsigma)
        MU <- -etapars[control$sigmaidx]^2/2
        gamma <- cholsigmainv%*%(Y-MU)  
        
        logpost <- logPosterior(surv=surv,X=X,beta=beta,omega=omega,eta=eta,gamma=gamma,priors=priors,cov.model=cov.model,u=u,control=control)        
                          
        return(logpost)
    }

    npts <- 20
    if(leneta>=3){
        npts <- 10
    }
    rgs <- getparranges(priors=priors,leneta=leneta)   
    qa <- QuadApprox(logpost,npts=npts,argRanges=rgs,plot=control$plotcal,surv=surv,X=X,beta=betahat,omega=omegahat,Y=Yhat,priors=priors,cov.model=cov.model,u=u,control=control)
    
    matr <- qa$curvature
    etahat <- qa$max
    
    # entry for eta in proposal covariance
    sigma[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta),(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)] <- matr    
        
    #estimate of gamma
    etahatpars <- cov.model$itrans(etahat)  
    ssigma <- matrix(EvalCov(cov.model=cov.model,u=u,parameters=etahatpars),n,n)
    cholssigma <- t(chol(ssigma))
    MU <- -etahatpars[control$sigmaidx]^2/2
    gammahat <- solve(cholssigma)%*%(Yhat-MU)
    
    hessian <- logPosterior(surv=surv,X=X,beta=betahat,omega=omegahat,eta=etahat,gamma=gammahat,priors=priors,cov.model=cov.model,u=u,control=control,hessian=TRUE)
    
    # beta and omega
    sigma[1:lenbeta,1:lenbeta] <- hessian$hess_beta
    sigma[(lenbeta+1):(lenbeta+lenomega),(lenbeta+1):(lenbeta+lenomega)] <- hessian$hess_omega
    sigma[(lenbeta+1):(lenbeta+lenomega),(1:lenbeta)] <- hessian$hess_omega_beta
    sigma[(1:lenbeta),(lenbeta+1):(lenbeta+lenomega)] <- t(hessian$hess_omega_beta)       
    # gamma
    diag(sigma)[(lenbeta+lenomega+leneta+1):npars] <- hessian$hess_gamma   
    
    return(list(etahat=etahat,sigma=solve(-sigma))) 
}




##' proposalVariance_gridded function
##'
##' A function to compute an approximate scaling matrix for the MCMC algorithm. Not intended for general use. 
##'
##' @param X the design matrix, containing covariate information 
##' @param surv an object of class Surv 
##' @param betahat an estimate of beta 
##' @param omegahat an estimate of omega 
##' @param Yhat an estimate of Y 
##' @param priors the priors 
##' @param cov.model the spatial covariance model 
##' @param u a vector of pairwise distances 
##' @param control a list containg various control parameters for the MCMC and post-processing routines 
##' @return an estimate of eta and also an approximate scaling matrix for the MCMC
##' @export

proposalVariance_gridded <- function(X,surv,betahat,omegahat,Yhat,priors,cov.model,u,control){

    Ygrid <- gridY(Y=Yhat,control=control)    
     
    n <- nrow(X)
    lenbeta <- length(betahat)
    lenomega <- length(omegahat)
    leneta <- getleneta(cov.model)
    lenY <- length(Ygrid)
    npars <- lenbeta + lenomega + leneta + lenY
    
    # eta
    logpost <- function(eta,surv,X,beta,omega,Ygrid,priors,cov.model,u,control){
        
        etapars <- cov.model$itrans(eta)
        covbase <- matrix(EvalCov(cov.model=cov.model,u=u,parameters=etapars),control$Mext,control$Next)
        
        rootQeigs <- sqrt(1/Re(fft(covbase)))   
       
        ymean <- -etapars[control$sigmaidx]^2/2
        gamma <- GammafromY(Ygrid,rootQeigs=rootQeigs,mu=ymean)  
        
        logpost <- logPosterior_gridded(surv=surv,X=X,beta=beta,omega=omega,eta=eta,gamma=gamma,priors=priors,cov.model=cov.model,u=u,control=control)        
                          
        return(logpost)
    }

    npts <- 20
    if(leneta>=3){
        npts <- 10
    }
    rgs <- getparranges(priors=priors,leneta=leneta)   
    qa <- QuadApprox(logpost,npts=npts,argRanges=rgs,plot=control$plotcal,surv=surv,X=X,beta=betahat,omega=omegahat,Ygrid=Ygrid,priors=priors,cov.model=cov.model,u=u,control=control)
    
    matr <- qa$curvature
    etahat <- qa$max
    
    sigma <- matrix(0,lenbeta + lenomega + leneta,lenbeta + lenomega + leneta)
    
    # entry for eta in proposal covariance
    sigma[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta),(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)] <- matr    
        
    #estimate of gamma
    etahatpars <- cov.model$itrans(etahat)  
    covbase <- matrix(EvalCov(cov.model=cov.model,u=u,parameters=etahatpars),control$Mext,control$Next)        
    rootQeigs <- sqrt(1/Re(fft(covbase)))   
    ymean <- -etahatpars[control$sigmaidx]^2/2
    gammahat <- GammafromY(Ygrid,rootQeigs=rootQeigs,mu=ymean)
    
    hessian <- logPosterior_gridded(surv=surv,X=X,beta=betahat,omega=omegahat,eta=etahat,gamma=gammahat,priors=priors,cov.model=cov.model,u=u,control=control,hessian=TRUE)
    
    # beta and omega
    sigma[1:lenbeta,1:lenbeta] <- hessian$hess_beta
    sigma[(lenbeta+1):(lenbeta+lenomega),(lenbeta+1):(lenbeta+lenomega)] <- hessian$hess_omega
    sigma[(lenbeta+1):(lenbeta+lenomega),(1:lenbeta)] <- hessian$hess_omega_beta
    sigma[(1:lenbeta),(lenbeta+1):(lenbeta+lenomega)] <- t(hessian$hess_omega_beta)       
    # gamma
    hess_gam <- hessian$hess_gamma 
    
    sigma <- (-1) * sigma # variance is inverse of observed information    
    
    matidx <- (lenbeta+lenomega+leneta+1):npars
    matidx <- matrix(matidx,nrow=length(matidx),ncol=2) 

    sigmaret <- Matrix(0,npars,npars)
    sigmaret[1:(lenbeta+lenomega+leneta),1:(lenbeta+lenomega+leneta)] <- solve(sigma)
    sigmaret[matidx] <- -1/hess_gam   
    
    return(list(etahat=etahat,sigma=sigmaret)) 
}

##' proposalVariance_polygonal function
##'
##' A function to compute an approximate scaling matrix for the MCMC algorithm. Not intended for general use. 
##'
##' @param X the design matrix, containing covariate information 
##' @param surv an object of class Surv 
##' @param betahat an estimate of beta 
##' @param omegahat an estimate of omega 
##' @param Yhat an estimate of Y 
##' @param priors the priors 
##' @param cov.model the spatial covariance model 
##' @param u a vector of pairwise distances 
##' @param control a list containg various control parameters for the MCMC and post-processing routines 
##' @return an estimate of eta and also an approximate scaling matrix for the MCMC
##' @export

proposalVariance_polygonal <- function(X,surv,betahat,omegahat,Yhat,priors,cov.model,u,control){

    Ygrid <- gridY_polygonal(Y=Yhat,control=control)    
     
    lenbeta <- length(betahat)
    lenomega <- length(omegahat)
    leneta <- getleneta(cov.model)
    lenY <- length(Ygrid)
    npars <- lenbeta + lenomega + leneta + lenY

    sigma <- matrix(0,npars,npars)

    
    # eta
    logpost <- function(eta,surv,X,beta,omega,Ygrid,priors,cov.model,u,control){

        etapars <- cov.model$itrans(eta)
        sigma <- matrix(EvalCov(cov.model=cov.model,u=u,parameters=etapars),control$n,control$n)
        cholsigma <- t(chol(sigma))
        cholsigmainv <- solve(cholsigma)
        MU <- -etapars[control$sigmaidx]^2/2
        gamma <- cholsigmainv%*%(Ygrid-MU)  
        
        logpost <- logPosterior_polygonal(surv=surv,X=X,beta=beta,omega=omega,eta=eta,gamma=gamma,priors=priors,cov.model=cov.model,u=u,control=control)        
                          
        return(logpost)
    }

    npts <- 20
    if(leneta>=3){
        npts <- 10
    }
    rgs <- getparranges(priors=priors,leneta=leneta)   
    qa <- QuadApprox(logpost,npts=npts,argRanges=rgs,plot=control$plotcal,surv=surv,X=X,beta=betahat,omega=omegahat,Ygrid=Ygrid,priors=priors,cov.model=cov.model,u=u,control=control)
    
    matr <- qa$curvature
    etahat <- qa$max
    
    # entry for eta in proposal covariance
    sigma[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta),(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)] <- matr    
        
    #estimate of gamma
    etahatpars <- cov.model$itrans(etahat)  
    ssigma <- matrix(EvalCov(cov.model=cov.model,u=u,parameters=etahatpars),control$n,control$n)
    cholssigma <- t(chol(ssigma))
    MU <- -etahatpars[control$sigmaidx]^2/2
    gammahat <- solve(cholssigma)%*%(Ygrid-MU)
    
    hessian <- logPosterior_polygonal(surv=surv,X=X,beta=betahat,omega=omegahat,eta=etahat,gamma=gammahat,priors=priors,cov.model=cov.model,u=u,control=control,hessian=TRUE)

    # if(any(eigen(hessian$hess_beta)$values<0)){
    #     cat("Problem with calibrating beta component, fixing using ad-hoc method...\n")
    #     hessian$hess_beta <- diag(abs(hessian$hess_beta))
    #     hessian$hess_omega_beta <- 0
    # }
    # if(any(eigen(hessian$hess_omega)$values<0)){
    #     cat("Problem with calibrating omega component, fixing using ad-hoc method ...\n")
    #     hessian$hess_omega <- diag(abs(hessian$hess_omega))
    #     hessian$hess_omega_beta <- 0
    # }

    # beta and omega
    sigma[1:lenbeta,1:lenbeta] <- hessian$hess_beta
    sigma[(lenbeta+1):(lenbeta+lenomega),(lenbeta+1):(lenbeta+lenomega)] <- hessian$hess_omega
    sigma[(lenbeta+1):(lenbeta+lenomega),(1:lenbeta)] <- hessian$hess_omega_beta
    sigma[(1:lenbeta),(lenbeta+1):(lenbeta+lenomega)] <- t(hessian$hess_omega_beta)       
    # gamma
    diag(sigma)[(lenbeta+lenomega+leneta+1):npars] <- hessian$hess_gamma   
    
    #browser()

    return(list(etahat=etahat,sigma=solve(-sigma))) 
}



##' proposalVariance_SPDE function
##'
##' A function to compute an approximate scaling matrix for the MCMC algorithm. Not intended for general use. 
##'
##' @param X the design matrix, containing covariate information 
##' @param surv an object of class Surv 
##' @param betahat an estimate of beta 
##' @param omegahat an estimate of omega 
##' @param Yhat an estimate of Y 
##' @param priors the priors 
##' @param cov.model the spatial covariance model 
##' @param u a vector of pairwise distances 
##' @param control a list containg various control parameters for the MCMC and post-processing routines 
##' @return an estimate of eta and also an approximate scaling matrix for the MCMC
##' @export

proposalVariance_SPDE <- function(X,surv,betahat,omegahat,Yhat,priors,cov.model,u,control){

    Ygrid <- gridY_polygonal(Y=Yhat,control=control)    


     
    lenbeta <- length(betahat)
    lenomega <- length(omegahat)
    leneta <- getleneta(cov.model)
    lenY <- length(Ygrid)
    npars <- lenbeta + lenomega + leneta + lenY

    sigma <- matrix(0,npars,npars)
    
    # eta
    logpost <- function(eta,surv,X,beta,omega,Ygrid,priors,cov.model,u,control){

        etapars <- cov.model$itrans(eta)
        prec <- (1/etapars[1])*control$precmat(SPDEprec(etapars[2],cov.model$order))
        U <- Matrix::chol(prec)

        if(cov.model$order>1){
            margvar <- etapars[1]/(4*pi*(cov.model$order-1)*(etapars[2]-4)^(cov.model$order-1)) # marginal var of Y = psi*variance_of_GMRF
        }
        else{
            margvar <- etapars[1]*etapars[2]/(4*pi)
        }
        MU <- rep(-margvar/2,control$n)
        gamma <- GammaFromY_SPDE(Ygrid,U=U,mu=MU)

        # y <- YFromGamma_SPDE(gamma=gamma,U=U,mu=MU)
        # hist(Ygrid-y)
        # browser()

        logpost <- logPosterior_SPDE(surv=surv,X=X,beta=beta,omega=omega,eta=eta,gamma=gamma,priors=priors,cov.model=cov.model,u=u,control=control)        
                          
        return(logpost)
    }

    npts <- 20
    if(leneta>=3){
        npts <- 10
    }
    rgs <- getparranges(priors=priors,leneta=leneta)   
    qa <- QuadApprox(logpost,npts=npts,argRanges=rgs,plot=control$plotcal,surv=surv,X=X,beta=betahat,omega=omegahat,Ygrid=Ygrid,priors=priors,cov.model=cov.model,u=u,control=control)
    
    matr <- qa$curvature
    etahat <- qa$max
    
    # entry for eta in proposal covariance
    sigma[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta),(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)] <- matr    
        
    #estimate of gamma
    etahatpars <- cov.model$itrans(etahat)  
    prec <- (1/etahatpars[1])*control$precmat(SPDEprec(etahatpars[2],cov.model$order))
    # cholprec <- t(chol(prec))
    # margvar <- etahatpars[1]/(4*pi*(cov.model$order-1)*(etahatpars[2]-4)^(cov.model$order-1)) # marginal var of Y = psi*variance_of_GMRF
    # MU <- rep(-margvar/2,control$n)
    # gammahat <- GammaFromY_SPDE(Ygrid,P=prec,L=cholprec,mu=MU)
    U <- chol(prec)
    if(cov.model$order>1){
        margvar <- etahatpars[1]/(4*pi*(cov.model$order-1)*(etahatpars[2]-4)^(cov.model$order-1)) # marginal var of Y = psi*variance_of_GMRF
    }
    else{
        margvar <- etahatpars[1]*etahatpars[2]/(4*pi)
    }
    MU <- rep(-margvar/2,control$n)
    gammahat <- GammaFromY_SPDE(Ygrid,U=U,mu=MU)
    
    hessian <- logPosterior_SPDE(surv=surv,X=X,beta=betahat,omega=omegahat,eta=etahat,gamma=gammahat,priors=priors,cov.model=cov.model,u=u,control=control,hessian=TRUE)
    
    # beta and omega
    sigma[1:lenbeta,1:lenbeta] <- hessian$hess_beta
    sigma[(lenbeta+1):(lenbeta+lenomega),(lenbeta+1):(lenbeta+lenomega)] <- hessian$hess_omega
    sigma[(lenbeta+1):(lenbeta+lenomega),(1:lenbeta)] <- hessian$hess_omega_beta
    sigma[(1:lenbeta),(lenbeta+1):(lenbeta+lenomega)] <- t(hessian$hess_omega_beta)       
    # gamma
    diag(sigma)[(lenbeta+lenomega+leneta+1):npars] <- hessian$hess_gamma   
    
    return(list(etahat=etahat,sigma=solve(-sigma))) 
}



##' estimateY function
##'
##' A function to get an initial estimate of Y, to be used in calibrating the MCMC. Not for general use
##'
##' @param X the design matrix containing covariate information 
##' @param betahat an estimate of beta 
##' @param omegahat an estimate of omega 
##' @param surv an object of class Surv
##' @param control a list containg various control parameters for the MCMC and post-processing routines 
##' @return an estimate of Y, to be used in calibrating the MCMC
##' @export

estimateY <- function(X,betahat,omegahat,surv,control){
   
    omega <- control$omegaitrans(omegahat) # this is omega on the correct scale
    
    haz <- setupHazard(dist=control$dist,pars=omega,grad=FALSE,hess=FALSE)    
    
    tsubs <- guess_t(surv)   
    
    Y <- -X%*%betahat - log(haz$H(tsubs)) # greedy estimate of Y (maximise individual contributions to log-likelihood) ... note log(delta) is now omitted  

    return(Y)    
}



##' guess_t function
##'
##' A function to get an initial guess of the failure time t, to be used in calibrating the MCMC. Not for general use
##'
##' @param surv an object of class Surv 
##' @return a guess at the failure times
##' @export

guess_t <- function(surv){

    n <- nrow(surv)
    
    censoringtype <- attr(surv,"type")
    
    if(censoringtype=="left" | censoringtype=="right"){
        notcensored <- surv[,"status"]==1
    }
    else{
        rightcensored <- surv[,"status"] == 0
        notcensored <- surv[,"status"] == 1
        leftcensored <- surv[,"status"] == 2
        intervalcensored <- surv[,"status"] == 3
    }   

    # setup function J=exp(X%*%beta + Y)*H_0(t)
    if(censoringtype=="left" | censoringtype=="right"){
        tsubs <- surv[,"time"]
    }
    else{ # else interval censored
        tsubs <- surv[,"time1"]        
    } 
    
    for(i in 1:n){
    
        if(notcensored[i]){
            next
        }

        if(censoringtype=="left" | censoringtype=="right"){            
            if(censoringtype=="right"){
                tpot <- tsubs[notcensored][tsubs[notcensored]>tsubs[i]] # potential t                    
            }
            else{ # censoringtype=="left"
                tpot <- tsubs[notcensored][tsubs[notcensored]<tsubs[i]] # potential t
            }
        }
        else{
            if(rightcensored[i]){
                tpot <- tsubs[notcensored][tsubs[notcensored]>tsubs[i]] # potential t
            }
            if(leftcensored[i]){
                tpot <- tsubs[notcensored][tsubs[notcensored]<tsubs[i]] # potential t
            }
            if(intervalcensored[i]){
                tpot <- surv[,"time1"][i] + 0.5*(surv[,"time2"][i]-surv[,"time1"][i]) # mid point of interval
            }
        }
            
        if(length(tpot)==0){
            next # leave tsubs[i] alone 
        }
        else{
            if(length(tpot)==1){
                tpot <- c(tpot,tpot)
            }
            tsubs[i] <- sample(tpot,1) # ignoring covariates, sample from empirical distribution of times exceeding (right censored), or less than (left censored) the observed time 
        }

    }
    
    return(tsubs)

}