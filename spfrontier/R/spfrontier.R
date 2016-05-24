# Copyright 2014 by Dmitry Pavlyuk <Dmitry.V.Pavlyuk@gmail.com>

#
# Spatial stochastic frontier model estimation
#


#' @title Spatial stochastic frontier model
#'
#' @description
#' \code{spfrontier} estimates spatial specifications of 
#' the stochastic frontier model.
#' 
#' @details
#' Models for estimation are specified symbolically, but without any spatial components.
#' Spatial components are included implicitly on the base of the \code{model} argument.
#' 
#' 
#' @param formula an object of class "\code{\link{formula}}": 
#' a symbolic description of the model to be fitted. 
#' The details of model specification are given under 'Details'.
#' @param data data frame, containing the variables in the model
#' @param W_y a spatial weight matrix for spatial lag of the dependent variable
#' @param W_v a spatial weight matrix for spatial lag of the symmetric error term
#' @param W_u a spatial weight matrix for spatial lag of the inefficiency error term
#' @param initialValues an optional vector of initial values, used by maximum likelihood estimator.
#' If not defined, estimator-specific method of initial values estimation is used.
#' @param logging an optional level of logging. Possible values are 'quiet','warn','info','debug'. 
#' By default set to quiet.
#' @param inefficiency sets the distribution for inefficiency error component. Possible values are 'half-normal' (for half-normal distribution) and 'truncated' (for truncated normal distribution). 
#' By default set to 'half-normal'. See references for explanations
#' @param onlyCoef allows calculating only estimates for coefficients (with inefficiencies and other additional statistics). Developed generally for testing, to speed up the process.
#' @param costFrontier is designed for selection of cost or production frontier
#' @param control an optional list of control parameters, 
#' passed to \code{\link{optim}} estimator from the '\code{\link{stats}} package
#' 
#' 
#' @keywords spatial stochastic frontier
#' @export
#' @references 
#' Kumbhakar, S.C. and Lovell, C.A.K (2000), Stochastic Frontier Analysis, Cambridge University Press, U.K.
#' @examples
#' 
#' data( airports )
#' airports2011 <- subset(airports, Year==2011)
#' W <- constructW(cbind(airports2011$longitude, airports2011$latitude),airports2011$ICAO)
#' formula <- log(PAX) ~ log(Population100km) + log(Routes) + log(GDPpc)
#
#' ols <- lm(formula , data=airports2011)
#' summary(ols )
#' plot(density(stats::residuals(ols)))
#' skewness(stats::residuals(ols))
#' 
#' # Takes >5 sec, see demo for more examples
#' # model <- spfrontier(formula , data=airports2011)
#' # summary(model )
#' 
#' # model <- spfrontier(formula , data=airports2011, W_y=W)
#' # summary(model )
#' 
spfrontier <- function(formula, data,
                       W_y = NULL, W_v = NULL,W_u = NULL,
                       inefficiency = "half-normal",
                       initialValues="errorsarlm",
                       logging = c("quiet", "info", "debug"),
                       control=NULL,
                       onlyCoef = F,
                       costFrontier = F){
    #Validation of model parameters
    start <- Sys.time()
    logging <- match.arg(logging)
    if (is.null(initialValues)) initialValues<-"errorsarlm"
    con <- list(grid.beta0 = 1, grid.sigmaV = 1, grid.sigmaU = 1, grid.rhoY = 1, grid.rhoU = 7, grid.rhoV = 7, grid.mu = 1,
                optim.control = list(tol=1e-5, iterlim=2000,repeatNM=1))
    namc <- names(control)
    con[namc] <- control
    noNms <- namc[!namc %in% names(con)]
    if (length(noNms)>0) 
        warning("Unknown parameters in control: ", paste(noNms, collapse = ", "))
    
    #End of model parameters' validation 
    
    initEnvir(W_y=W_y, W_v=W_v,W_u=W_u,inefficiency=inefficiency,
              initialValues=initialValues,logging=logging,costFrontier=costFrontier,control=con)
    logging("Estimator started", level="info")
    
    #Preparing the environment
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
    isSpY <- !is.null(W_y)
    isSpV <- !is.null(W_v)
    isSpU <- !is.null(W_u)
    isTN <- (inefficiency == "truncated")
    costf <- ifelse(costFrontier, -1, 1)
    if ((isSpV || isSpU) &&(logging!="debug")){
        logging("MLE can take a long time for spatial lags in error components. This is recommended to use logging='debug' to control the progress", level="info")
    }
    #if(!is.null(tv)){
    #    logging("logL at true value:", funcLogL.direct(tv))
    #}
    
    if (typeof(initialValues)!="double"){
        logging("Calculating initial values",level="info")
        iniParams <- calculateInitialValues(formula, data)
        initialValues <- paramsToVector(olsenReparam(iniParams))
        logging("Initial values:", paramsToVector(iniParams, olsen=F),level="info")
    }else{
        initialValues <- paramsToVector(olsenReparam(paramsFromVector(initialValues, k, isSpY, isSpV, isSpU, isTN, olsen=F)))
    }
    grad <- NULL
    if (!isSpV && !isSpU){
        grad <- funcGradient
    }
    constr <- NULL
    if (isSpY || isSpV || isSpU){
        num <- 0
        if (isSpY) num <- num + 1
        if (isSpV) num <- num + 1
        if (isSpU) num <- num + 1
        vars <- k+num+2
        cons <- num*2 + 2
        if (isTN) vars <- vars + 1
        
        pscale <- rep(1, vars)
        A <- matrix(rep(0, vars*cons), ncol=vars,nrow=cons)
        B <- rep(0, cons)
        fr <- k
        if (isSpY) fr <- fr+1
        A[1,fr+1] <- 1
        A[2,fr+2] <- 1
        l <- 2
        if (isSpY){
            A[l+1,fr] <- -1
            A[l+2,fr]<- 1
            B[l+1] <- 1
            B[l+2] <- 1
            l <- l+2
        }
        if (isSpV){
            A[l+1,fr+3] <- -1
            A[l+2,fr+3]<- 1
            B[l+1] <- 1
            B[l+2] <- 1
            pscale[fr+3] <- n
            l <- l+2
            fr <- fr + 1
        }
        if (isSpU){
            A[l+1,fr+3] <- -1
            A[l+2,fr+3]<- 1
            pscale[fr+3] <- n
            B[l+1] <- 1
            B[l+2] <- 1
        }
        constr <- list(ineqA=A, ineqB=B)
        envirAssign("constr",constr)
        envirAssign("parscale",pscale)
    }
    estimates <- optimEstimator(formula, data, 
                                funcLogL, 
                                ini = initialValues, 
                                gr=grad)
    logging(paste("Completed, status =", status(estimates)), level = "info")
    if(status(estimates) == 0){
        olsenP <- paramsFromVector(resultParams(estimates), k, isSpY, isSpV, isSpU, isTN)
    p <- olsenReparamBack(olsenP)
    coefficients(estimates) <- p
    
    logging("Final estimates:", paramsToVector(p, olsen=F), level = "info")
    
    if(!onlyCoef){
        logging("Preparing additional statistics...")
        fittedY <- X %*% p$beta
        if (!is.null(p$rhoY)){
            Sp1 <- solve(diag(n) - p$rhoY * W_y)
            fittedY <- Sp1 %*% fittedY
        }
        rownames(fittedY) <- rownames(X)
        colnames(fittedY) <- c("Fitted values")
        fitted(estimates) <- fittedY
        
        resid <- y - fittedY
        rownames(resid) <- rownames(y)
        colnames(resid) <- c("Residuals")
        residuals(estimates) <- resid
    
        tryCatch({
            #logging("Calculating hessian")
            #hessian(estimates) <- numDeriv::hessian(funcLogL,x=resultParams(estimates))
            logging("Calculating stdErrors...")
            #hess <- optimHess(paramsToVector(p, olsen = F),funcLogL.direct)
            #hessian(estimates) <- hess
            
            G <- olsenGradient(olsenP)
            iH <- solve(-hessian(estimates))
            stdErrors(estimates) <- sqrt(diag(G%*%iH%*%t(G)))
            logging("Done")
        }, error = function(e){
            logging(e$message, level="warn")
        })
        tryCatch({
            logging("Calculating efficiencies...")
            mu <- 0
            if(isTN){
                mu <- p$mu
            }
            if (is.null(p$rhoV) && is.null(p$rhoU) ){
                sigma <- sqrt(p$sigmaU^2 + p$sigmaV^2)
                lambda <- p$sigmaU / p$sigmaV
                A <- costf*resid * lambda / sigma - mu / (lambda * sigma)
                u <- (dnorm(A) / (1 - pnorm(A)) - A) * p$sigmaU * p$sigmaV / sigma
            }else{
                I <- diag(n)
                SpV <- I
                SpU <- I
                
                if (!is.null(p$rhoV))
                    SpV <- solve(I-p$rhoV*W_v)
                if (!is.null(p$rhoU))
                    SpU <- solve(I-p$rhoU*W_u)
                
                mSigmaV = p$sigmaV^2*SpV%*%t(SpV)
                mSigmaU = p$sigmaU^2*SpU%*%t(SpU)
                mSigma = mSigmaV + mSigmaU
                imSigma <- solve(mSigma)
                mDelta = mSigmaU%*%imSigma%*%mSigmaV 
                rownames(mDelta) <- colnames(mDelta)
                
                mA = mDelta %*% solve(mSigmaV)
                mD = -mDelta %*% solve(mSigmaU)
                
                vMu = rep(mu, n)
                mGamma <- -costf*mSigmaU%*%imSigma
                m <- vMu + mGamma %*% (resid + costf*vMu)
                u <- mtmvnorm(lower= rep(0, n),mean=as.vector(t(m)), sigma=mDelta, doComputeVariance=FALSE)$tmean
            }
            eff <- cbind(exp(-u))
            rownames(eff) <- rownames(y)
            colnames(eff) <- c("Efficiency values")
            efficiencies(estimates) <- eff
        }, error = function(e){
            logging(e$message,p$rhoV, level="warn")
        })
    }
    }
    logging(paste("Done! Elapsed time: ",Sys.time()-start))
    finalizeEnvir()
    
    return(estimates)
}

calculateInitialValues <- function (formula, data) {
    y <- envirGet("y")
    X <- envirGet("X")
    k <- ncol(X)
    W_y <- envirGet("W_y")
    W_v <- envirGet("W_v")
    W_u <- envirGet("W_u")
    costFrontier <- envirGet("costFrontier")
    inefficiency <- envirGet("inefficiency")
    initialValues <- envirGet("initialValues")
    logging(paste("Initial values method",initialValues),level="info")
    isSpY <- !is.null(W_y)
    isSpV <- !is.null(W_v)
    isSpU <- !is.null(W_u)
    isTN <- (inefficiency == "truncated")
    iniParams <- list()
    f <- formula
    if (isSpY){
        data$Wy <- W_y %*% y
        f <- update(formula,    ~ . + Wy)
    }
    e <- NULL
    found <- F
    if (isSpY){
        iniParams$rhoY <- 0
        inconsistentSpatial <- spfrontier(f, data, logging="quiet",onlyCoef=F,costFrontier=costFrontier)
        coef <- coefficients(inconsistentSpatial)
        if (length(coef)>0){
            iniParams$rhoY = tail(coef$beta, n=1)
            iniParams$beta = head(coef$beta, -1)
            iniParams$sigmaU <- coef$sigmaU
            iniParams$sigmaV <- coef$sigmaV
            e <- residuals(inconsistentSpatial)
            found <- T
        }
    }
    if(!found){
        sfa <- sfaInitialValues(f, data)
        if (isSpY){
            iniParams$rhoY = tail(sfa$beta, n=1)
            iniParams$beta = head(sfa$beta, -1)
        }else{
            iniParams$beta <- sfa$beta
        }
        e <- sfa$residuals
        iniParams$sigmaU <- sfa$sigmaU
        iniParams$sigmaV <- sfa$sigmaV
    }
    
    mu <- -mean(e) 
    
    if (is.null(initialValues) || initialValues == "nonspatial"){
        if (isSpV){
            iniParams$rhoV <- 0
        }
        if (isSpU){
            iniParams$rhoU <- 0
        }
    }
    if (!is.null(initialValues)&&initialValues == "errorsarlm"){
        if (isSpV){
            listw <- mat2listw(W_v)
            sem <- errorsarlm(f ,data=data, listw)
            logging(paste("SEM lambda for W_v",sem$lambda),level="debug")
            iniParams$rhoV <- sem$lambda
            #We <- W_v %*% e
            #sarResiduals <- lm(e ~ We - 1, data = data.frame(e, We))
            #iniParams$rhoV <- coef(sarResiduals)
            #mu <- -mean(stats::residuals(sarResiduals))
        }
        if (isSpU){
            listw <- mat2listw(W_u)
            sem <- errorsarlm(f ,data=data, listw)
            logging(paste("SEM lambda for W_u",sem$lambda),level="debug")
            iniParams$rhoU <- sem$lambda
        }
    }
    if(isTN){
        iniParams$mu <- mu
    }
    if (!is.null(initialValues)  && initialValues=="grid" &&(isSpU || isSpV)){
        iniParams <- gridSearch(iniParams)
    }
    initialValues <- iniParams
    return(initialValues)
}

gridSearch <- function(params){
    con <- envirGet("control")
    
    if (is.null(params$rhoV) && is.null(params$rhoU)){
        return(params)
    }
    logging("Grid search for initial values")
    gridVal <- list()
    c <-0 
    for(b in params$beta){
        if (c == 0){
            gridVal[[paste("beta",c, sep="")]] <- seq(b,b + 3*params$sigmaV,length.out = con$grid.beta0)
        }else{
            gridVal[[paste("beta",c, sep="")]] <- b
        }
        c <- c+1
    }
    if (!is.null(params$rhoY)){
        gridVal$rho = params$rhoY
    }
    gridVal$sigmaV = c(params$sigmaV,seq(0.5*params$sigmaV, 3*params$sigmaV, length.out = con$grid.sigmaV))
    #gridVal$sigmaU = params$sigmaU
    gridVal$sigmaU = c(params$sigmaU,seq(0.5*params$sigmaU, 3*params$sigmaU, length.out = con$grid.sigmaU))
    if (!is.null(params$rhoV)){
        gridVal$rhoV = seq(-0.3, 0.3, length.out = con$grid.rhoV)
    }
    if (!is.null(params$rhoU)){
        gridVal$rhoU = seq(-0.3, 0.3, length.out = con$grid.rhoU)
    }
    if (!is.null(params$mu)){
        gridVal$mu = seq(params$mu-3*params$sigmaV,params$mu + 3*params$sigmaV,length.out = con$grid.mu)
    }
    grid <- expand.grid(gridVal)
    
    ml <- apply(grid, 1, funcLogL.direct)
    opt <- min(ml)
    ret <- as.vector(as.matrix(grid[which.min(ml), ]))
    p <- paramsFromVector(ret,length(params$beta),!is.null(params$rhoY),!is.null(params$rhoV),!is.null(params$rhoU),!is.null(params$mu), olsen = F)
    if (opt == 1e16){
        if (!is.null(p$rhoV)) p$rhoV<-0
        if (!is.null(p$rhoU)) p$rhoU<-0
    }
    return(p)
}

sfaInitialValues <- function(formula, data){
    ols <- lm(formula, data=data)
    res_ols <- stats::resid(ols)
    m2 = sum(res_ols^2)/length(res_ols)
    m3 = sum(res_ols^3)/length(res_ols)
    sigmaU = (m3*sqrt(pi/2)/(1-4/pi))^(1/3)
    sigmaV = 0
    if (!is.nan(sigmaU) && (m2 - (1-2/pi)*sigmaU^2>0)){
        sigmaV = sqrt(m2 - (1-2/pi)*sigmaU^2)
    }
    if ((sigmaU<=0) || (sigmaV<=0)){
        sigmaU <- var(res_ols)*(1-2/pi)
        sigmaV <- var(res_ols)
    }
    beta <- as.vector(coef(ols))
    
    return(list(beta = beta,sigmaV=sigmaV,sigmaU=sigmaU,residuals = res_ols))
}
