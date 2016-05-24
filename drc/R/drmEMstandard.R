"drmEMstandard" <- 
function(dose, resp, multCurves, doseScaling = 1)
{

    ## Defining a helper function for calculating the variance-covariance matrix
#    vcFct <- function(beta0, beta, sigma2, len0)
#    {
#        vc <- (sigma2 / len0) * (beta %o% beta) / (beta0^4)
#        diag(vc) <- diag(vc) + sigma2 / (beta0^2)
#
#        return(vc)
#    }
    vcFct <- function(beta0, beta, len0)
    {
        vc <- (1 / len0) * (beta %o% beta) / (beta0^4)
        diag(vc) <- diag(vc) + (1 / (beta0^2))

        return(vc)
    }
    
    zeroDose <- dose < 1e-15  # hardcoded tolerance of 1e-15
    len0 <- sum(zeroDose)
    vcFct2 <- function(beta0, betaVec)
    {
#        len0 <- weightVec[1]  # in case len0 is a vector
    
        vc <- (1 / len0) * (betaVec %o% betaVec) / (beta0^4)
        diag(vc) <- diag(vc) + (1 / (beta0^2))
        
#        zeroDose <- dose < doseTol
#        print(vc[!zeroDose, zeroDose])
#        print((1 / len0) * (-betaVec / (beta0^3)))
#        print(vc[!zeroDose, zeroDose] + (1 / len0) * (-betaVec[!zeroDose] / (beta0^3)))
        
        vc[!zeroDose, zeroDose] <- vc[!zeroDose, zeroDose] + (1 / len0) * (-betaVec[!zeroDose] / (beta0^3)) 
        vc[zeroDose, !zeroDose] <- vc[zeroDose, !zeroDose] + (1 / len0) * (-betaVec[!zeroDose] / (beta0^3))
#        print(vc[zeroDose, zeroDose])
#        print(diag(vc[zeroDose, zeroDose]) + (1 / (len0 * beta0^2)) - (1 / (beta0^2)))
        diag(vc[zeroDose, zeroDose]) <- diag(vc[zeroDose, zeroDose]) + (1 / (len0 * beta0^2)) - (1 / (beta0^2))

        return(vc)
    }    


    ## Defining the objective function                
    opfct <- function(c)  # dose, resp and weights are fixed
    {
        print(c) 
        f0 <- multCurves(0, c)[1]
        print(f0)
        fVec <- multCurves(dose / doseScaling, c)
        print(fVec)
#        vcMat <- vcFct(f0, fVec, weightVec)   
        vcMat <- vcFct2(f0, fVec)   
        print(solve(vcMat)[1:6, 1:6])
        
        sum( (resp - fVec) %*% solve(vcMat) %*% (resp - fVec))
    }    

    
    ## Defining self starter function
    ssfct <- NULL


    ## Defining the log likelihood function
    llfct <- function(object)
    {
#        total <- (object$"data")[iv, 5]
#        success <- total*(object$"data")[iv, 2]    
#        c( sum(log(choose(total, success))) - object$"fit"$"ofvalue", object$"sumList"$"df.residual" )
        
        c(
        -object$"fit"$value + sum(log(gamma(resp+1))),
        object$"sumList"$"df.residual"
        )  # adding scale constant
    }
    
       
    ## Defining functions returning the residual variance, the variance-covariance matrix, and the parameter estimates
#    rvfct <- function(object)
#    {
#        object$"fit"$"value" / df.residual(object)  # object$"sumList"$"df.residual"
#    }
#
#    vcovfct <- function(object)
#    {
#        solve(object$fit$hessian)    
#    }
#

    # copied from drmEMls.R
    rvfct <- function(object)
    {
        object$"fit"$"value" / df.residual(object)
    }

    vcovfct <- function(object)
    {
        scaledH <- (object$"fit"$"hessian") / (2 * rvfct(object))
        invMat <- try(solve(scaledH), silent = TRUE)
    
        if (inherits(invMat, "try-error"))
        {
            ## More stable than 'solve' (suggested by Nicholas Lewin-Koh - 2007-02-12)
            ch <- try(chol(scaledH))
            if(inherits(ch, "try-error")) 
            {
                ch <- try(chol(0.99 * object$fit$hessian + 0.01 * diag(dim(object$fit$hessian)[1])))
            }
            ## Try regularizing if the varcov is unstable
            if(!inherits(ch, "try-error")) return(chol2inv(ch))
        } else {
            return(invMat)
        }
    } 
    
    parmfct <- function(fit, fixed = TRUE)
    {
        fit$par
    }


    ## Returning list of functions
    return(list(llfct = llfct, opfct = opfct, ssfct = ssfct, rvfct = rvfct, vcovfct = vcovfct, 
    parmfct = parmfct))
}


"drmLOFstandard" <- function()
{
    return(list(anovaTest = NULL, gofTest = NULL))
}



if (FALSE)
{


covFct <- function(sigma0, sigma, myVec, dose)
{
    zeroDose <- dose < 1e-15  # hardcoded tolerance of 1e-15
    len0 <- sum(zeroDose)
    
    my0 <- (myVec[zeroDose])[1]
    n0 <- sum(zeroDose)

    lenMy <- length(myVec)
    derMat <- matrix(0, lenMy, lenMy+1)
    diag(derMat) <- 1 / my0
    derMat[, lenMy+1] <- -myVec / (my0^2)

    lenY <- lenMy + 1
    origVCmat <- matrix(0, lenY, lenY)
    sigma0mean <- sigma0^2/n0
    origVCmat[zeroDose, zeroDose] <- sigma0mean
    
    diag(origVCmat)[!zeroDose] <- sigma^2
    diag(origVCmat)[zeroDose] <- sigma0^2
    
    origVCmat[lenY, lenY] <- sigma0mean
       
    list(my0, n0, derMat, origVCmat, derMat %*% origVCmat %*% t(derMat))
}


resList<-covFct(0.52, 0.52, fitted(ryegrass.m1)[1:10], ryegrass$conc[1:10])
resList[[5]] / (outVec %o% outVec)


varOptim1 <- function(varpar)
{
    resList <- covFct(varpar[1], varpar[2], fitted(ryegrass.m1), ryegrass$conc)[[5]]
    resVec <- residuals(ryegrass.m1)
#    resVec <- fitted(ryegrass.m1) + rnorm(24, 0, c(rep(3,6), rep(1, 18))
    resVec%*%solve(resList)%*%resVec + log(abs(det(resList)))

}

optim(c(1, 0.1), varOptim1)


varOptim1b <- function(par, const = 1)
{
    fittedVec <- par[2]+(1-par[2])/(1+(ryegrass$conc/par[3])^par[1])
    resList <- covFct(1, par[4], fittedVec, ryegrass$conc)[[5]]
    resVec <- (ryegrass$rootl / 7.75) - fittedVec
    resVec%*%solve(resList)%*%resVec + const * log(abs(det(resList)))
}

rg.optim <- optim(c(2,0.05,3,0.5), varOptim1b, hessian = TRUE)

sqrt(diag(solve(rg.optim$hessian)))

sqrt(varOptim1b(rg.optim$par, 0) / 20)



## S.alba
varOptim1b2 <- function(par, const = 1)
{
    fittedVec <- par[2]+(1-par[2])/(1+(S.alba$Dose[1:32]/par[3])^par[1])
    resList <- covFct(1, par[4], fittedVec, S.alba$Dose[1:32])[[5]]
    resVec <- (S.alba$DryMatter[1:32] / 7.75) - fittedVec
    resVec%*%solve(resList)%*%resVec + const * log(abs(det(resList)))
}

rg.optim2 <- optim(c(2,0.05,3,0.5), varOptim1b2, hessian = TRUE)
rg.optim2$par

sqrt(diag(solve(rg.optim2$hessian)))

sqrt(varOptim1b2(rg.optim2$par, 0) / 28)

sa.drm1 <- drm(DryMatter~Dose, data=S.alba[1:32,], fct=LL.4())
summary(sa.drm1)
sa.drm2 <- drm(DryMatter/mean(S.alba$DryMatter[1:8])~Dose, data=S.alba[1:32,], fct=LL.4(fixed=c(NA,NA,1,NA)))
summary(sa.drm2)
sa.drm3 <- drm(DryMatter/mean(S.alba$DryMatter[1:8])~Dose, data=S.alba[1:32,], fct=LL.4(fixed=c(NA,NA,NA,NA)))
summary(sa.drm3)


#yVec <- fitted(ryegrass.m1) + rnorm(24, 0, c(rep(3,6), rep(1, 18)))
yVec <- fitted(ryegrass.m1) + rnorm(24, 0, c(rep(10,6), rep(1, 18)))
xVec <- ryegrass$conc
varOptim1c <- function(par, const = 1)
{
    fittedVec <- par[2]+(1-par[2])/(1+(xVec/par[3])^par[1])
    resList <- covFct(1, par[4], fittedVec, xVec)[[5]]
    resVec <- (yVec /  mean(yVec[1:4])) - fittedVec
    resVec%*%solve(resList)%*%resVec + const * log(det(resList))
}

ratioVec <- rep(NA, 100)
sigmaVec <- rep(NA, 100)
ec50Vec <- rep(NA, 100)
seVec1 <- rep(NA, 100)
seVec2 <- rep(NA, 100)

xVec <- rep(ryegrass$conc, rep(3, 24))
xVec <- xVec[-c(1:14)]
for (i in 1:100)
{
    yVec <- rep(fitted(ryegrass.m1), rep(3, 24)) + rnorm(72, 0, c(rep(3,18), rep(1, 54)))
    yVec <- yVec[-c(1:14)]

#    varOptim1c <- function(par, const = 1)
#    {
#        fittedVec <- par[2]+(1-par[2])/(1+(xVec/par[3])^par[1])
#        resList <- covFct(1, par[4], fittedVec, xVec)[[5]]
#        resVec <- (yVec /  mean(yVec[1:6])) - fittedVec
#        resVec%*%solve(resList)%*%resVec + const * log(det(resList))
#    }  
#    seVec1[i] <- coef(summary(drm(yVec/mean(yVec[1:4])~xVec, fct=LL.4(fixed=c(NA,NA,NA,NA)))))[4,2]
    seVec1[i] <- coef(summary(drm(yVec/mean(yVec[1:4])~xVec, fct=LL.4(fixed=c(NA,NA,1,NA)))))[3,2]
    optimRes <- optim(c(1.7,0.01,2.5,0.1), varOptim1c, hessian=TRUE)
    seVec2[i] <- sqrt(diag(solve(optimRes$hessian)))[3]
    parem <- optimRes$par    
    ratioVec[i] <- parem[4]
    ec50Vec[i] <- parem[3]
    sigmaVec[i] <- sqrt(varOptim1c(parem, 0)/68)
}

cbind(seVec1, seVec2)

ratioVec
sigmaVec
hist(ratioVec)
hist(sigmaVec)

    




varOptim2 <- function(varpar)
{
    resList<-covFct(varpar[1], varpar[2], fitted(ryegrass.m1), ryegrass$conc)[[5]]
    resVec <-residuals(ryegrass.m1)
    resVec%*%solve(resList)%*%resVec+log(abs(det(resList)))

}

#resVec2 <- ryegrass$rootl - (fitted(ryegrass.m1) + rnorm(24, 0, c(rep(3,6), rep(1, 18))))
varOptim3 <- function(varpar, const=1)
{
    resList<-covFct(1, varpar[1], fitted(ryegrass.m1), ryegrass$conc)[[5]]
    resVec <- residuals(ryegrass.m1)
#    resVec <- resVec2
    resVec%*%solve(resList)%*%resVec + const * log(abs(det(resList)))

}

optimize(varOptim3, lower=0, upper=100)



}