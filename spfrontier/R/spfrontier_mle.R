Infin <- -1e16
funcLogL <- function(parameters, quiet = F){
    y <- envirGet("y")
    X <- envirGet("X")
    k <- ncol(X)
    W_y <- envirGet("W_y")
    W_v <- envirGet("W_v")
    W_u <- envirGet("W_u")
    costFrontier <- envirGet("costFrontier")
    inefficiency <- envirGet("inefficiency")
    
    isSpY <- !is.null(W_y)
    isSpV <- !is.null(W_v)
    isSpU <- !is.null(W_u)
    isTN <- (inefficiency == "truncated")
    costf <- ifelse(costFrontier, -1, 1)
          
    p <- paramsFromVector(parameters, k, isSpY, isSpV, isSpU, isTN)
    if(!quiet){
        pbeta <- olsenReparamBack(p)
        paramsToVector(pbeta, olsen = F)
        counter = envirCounter("funcLogL")
        logging(paste("Evaluating 'funcLogL' for parameters","[run=",counter,"]:"),paramsToVector(pbeta, olsen = F)) 
    }  
    n <- length(y)
    I <- diag(n)
    ret <- -Inf
    if ((p$nu>0) && (p$lambda>0) && (is.null(p$rhoY) || abs(p$rhoY)<1) 
        && (is.null(p$rhoV) || abs(p$rhoV)<1) 
        && (is.null(p$rhoU) || abs(p$rhoU)<1)){
        
        e <- y - X %*% p$gamma/p$nu
        
        if (isSpY){
            e <- e - p$rhoY * W_y %*% y
        }
        tryCatch({
            sigmaV <- 1/(p$nu*sqrt(1+p$lambda^2))
            sigmaU <- p$lambda * sigmaV
            mu <- 0
            if(isTN){
                mu <- p$mu
            }
            vMu = rep(mu, n)
            r <- 0
            if (isSpV || isSpU){
                SpV <- I
                SpU <- I
                if (isSpV)
                    SpV <- solve(I-p$rhoV*W_v)
                if (isSpU)
                    SpU <- solve(I-p$rhoU*W_u)
                mSigmaV = sigmaV^2*SpV%*%t(SpV)
                mSigmaU = sigmaU^2*SpU%*%t(SpU)
                mSigma = mSigmaV + mSigmaU
                imSigma <- solve(mSigma)
                mDelta <- mSigmaU%*%imSigma%*%mSigmaV 
                mGamma <- -costf*mSigmaU%*%imSigma
                rownames(mDelta) <- colnames(mDelta)
                
                f1 <- pmvnorm(lower=-Inf,upper = 0, mean=-vMu, sigma=mSigmaU)
                r <- -log(f1)
                r <- r + dmvnorm(x=as.vector(e),mean=-costf*vMu, sigma=mSigma,log=TRUE)
                f2 <- pmvnorm(lower=-Inf,upper = as.vector(mGamma%*%(e+costf*vMu)), mean=-vMu, sigma=mDelta)
                r <- r + log(f2)
            }else{
                r <- r - n*log(2*pi)/2
                r <- r - n*log(pnorm(mu/sigmaU))
                r <- r + n * log(p$nu) 
                r <- r - 0.5 * p$nu^2*sum((e + costf*vMu)^2)
                r <- r + sum(log(pnorm(mu*p$nu/p$lambda - costf*p$lambda*e*p$nu)))
            }
            if (isSpY)
                r <- r + determinant(I-p$rhoY*W_y,logarithm = TRUE)$modulus
            ret <- r
        }, error = function(e){
            logging(e$message, level="debug")
        })
    }else{
        if(!quiet) logging("Parameters are out of space")
    }
    if(!quiet) logging(paste("LogL =",ret,"\n"))
    #if (ret == -Inf || is.nan(ret)) ret<- Infin  For fake finite differences
    return(as.numeric(ret))
}

funcLogL.direct <- function(parameters){
    y <- envirGet("y")
    X <- envirGet("X")
    k <- ncol(X)
    W_y <- envirGet("W_y")
    W_v <- envirGet("W_v")
    W_u <- envirGet("W_u")
    inefficiency <- envirGet("inefficiency")
    
    isSpY <- !is.null(W_y)
    isSpV <- !is.null(W_v)
    isSpU <- !is.null(W_u)
    isTN <- (inefficiency == "truncated")
    
    p <- paramsFromVector(parameters, k, isSpY, isSpV, isSpU, isTN, olsen = F)
    polsen <- olsenReparam(p)
    return(funcLogL(paramsToVector(polsen), quiet = F))
}

funcGradient<-function(parameters){
    logging("Evaluating 'funcGradient' for parameters:",parameters)
    y = envirGet("y")
    X = envirGet("X")
    k <- ncol(X)
    W_y <- envirGet("W_y")
    W_v <- envirGet("W_v")
    W_u <- envirGet("W_u")
    inefficiency <- envirGet("inefficiency")
    costFrontier <- envirGet("costFrontier")
    
    isSpY <- !is.null(W_y)
    isSpV <- !is.null(W_v)
    isSpU <- !is.null(W_u)
    isTN <- (inefficiency == "truncated")
    costf <- ifelse(costFrontier, -1, 1)
    
    p <- paramsFromVector(parameters, k, isSpY, isSpV, isSpU, isTN)
    
    yst <- y
    if(isSpY){
        yst <- y-p$rhoY*W_y%*%y
    }
    mu <- 0
    if (isTN){
        mu <- p$mu
    }
    omega <- p$nu * yst - X%*%p$gamma
    a <- -costf*omega * p$lambda + mu*p$nu/p$lambda
    delta <- dnorm(a)/pnorm(a)
    N <- length(y)
    
    sigmaU <- p$lambda/(p$nu*sqrt(1+p$lambda^2))
    a2 <- mu/sigmaU
    delta2 <- (dnorm(a2)/pnorm(a2))
    
    dLdGamma <- t(omega+costf*mu*p$nu)%*%X + costf*t(delta)%*%X * p$lambda
    dLdNu <- -t(omega+costf*mu*p$nu)%*%(y+mu) - t(delta)%*%(costf*y*p$lambda-mu/p$lambda) + N * 1/p$nu - N*delta2*a2/p$nu
    dLdLambda <- -t(delta)%*%(costf*omega+mu*p$nu/p$lambda^2)+N*delta2*mu*p$nu/(p$lambda^2*sqrt(1+p$lambda^2))
    grad = c(dLdGamma)
    names(grad) = paste("dGamma", seq(length(dLdGamma)), sep = "")
    if(isSpY){
        mat = solve(diag(N)-p$rhoY*W_y) %*% (-W_y)
        dLdRhoY = p$nu*(t(omega)+costf*p$lambda*t(delta))%*%W_y%*%y+sum(diag(mat))
        ns <- names(grad)
        grad <- c(grad, dLdRhoY)
        names(grad) <- c(ns, "dRhoY")
    }
    ns <- names(grad)
    grad = c(grad, dLdNu,dLdLambda)
    names(grad) = c(ns, "dNu", "dLambda")
    if(isTN){
        dLdMu <- -N*delta2/sigmaU - p$nu*sum(omega+costf*mu*p$nu)+sum(delta)*p$nu/p$lambda
        #dLdMu = -(mu - 1)
        ns <- names(grad)
        grad <- c(grad, dLdMu)
        names(grad) <- c(ns, "dMu")
    }
    logging("Gradient",grad)
    return(grad)
}

#' @title Calculation of the log likelihood function for the spatial stochastic frontier model
#'
#' @description
#' \code{logLikelihood} returns a value of the log likelihood function 
#' for the spatial stochastic frontier model
#' 
#' @details
#' This function is exported from the package for testing and presentation purposes
#' A list of arguments of the function exactly matches the corresponding list of the \code{\link{spfrontier}} function
#' 
#' 
#' @param formula an object of class "\code{\link{formula}}"
#' @param data data frame, containing the variables in the model
#' @param W_y a spatial weight matrix for spatial lag of the dependent variable
#' @param W_v a spatial weight matrix for spatial lag of the symmetric error term
#' @param W_u a spatial weight matrix for spatial lag of the inefficiency error term
#' @param values a  vector of log likelihood function parameters
#' @param logging an optional level of logging. Possible values are 'quiet','warn','info','debug'. 
#' By default set to quiet.
#' @param inefficiency sets the distribution for inefficiency error component. Possible values are 'half-normal' (for half-normal distribution) and 'truncated' (for truncated normal distribution). 
#' By default set to 'half-normal'.
#' @param costFrontier is designed for selection of cost or production frontier

#' 
#' @export
logLikelihood <- function(formula, data,
                       W_y = NULL, W_v = NULL,W_u = NULL,
                       inefficiency = "half-normal",
                       values,
                       logging = c("quiet", "info", "debug"), costFrontier=F){
    logging <- match.arg(logging)
    con <- list(grid.beta0 = 1, grid.sigmaV = 1, grid.sigmaU = 1, grid.rhoY = 1, grid.rhoU = 10, grid.rhoV = 10, grid.mu = 1)
    
    initEnvir(W_y=W_y, W_v=W_v,W_u=W_u,inefficiency=inefficiency, logging=logging)
    logging("Estimator started", level="info")

    mf <- model.frame(formula, data)
    y <- as.matrix(model.response(mf))
    X <- as.matrix(mf[-1])
    tm <- attr(mf, "terms")
    intercept <- attr(tm, "intercept") == 1
    if (intercept){
        X <- cbind(Intercept=1L,X)
    }
    k <- ncol(X)
    n <- length(y)
    envirAssign("X", X)
    envirAssign("y", y)
    envirAssign("costFrontier", costFrontier)
    isSpY <- !is.null(W_y)
    isSpV <- !is.null(W_v)
    isSpU <- !is.null(W_u)
    isTN <- (inefficiency == "truncated")
    res <- funcLogL(paramsToVector(olsenReparam(paramsFromVector(values, k, isSpY, isSpV, isSpU, isTN, olsen=F))))
    names(res)<-"Log-likelihood function"
    return(res)
}