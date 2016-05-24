"drmEMls" <- 
function(dose, resp, multCurves, startVec, robustFct, weights, rmNA, dmf = NULL, 
doseScaling = 1, respScaling = 1)
{
#    ## Defining lack-of-fit/goodness-of-fit tests
#    anovaTest <- contAnovaTest()
#    gofTest <- NULL
#    if (anovaYes) {return(list(anovaTest = anovaTest, gofTest = gofTest))}

    ## Defining the objective function and its derivative
    opfct <- function(parm)  # , scaling = TRUE)
    {
#        print(parm)

#        if (scaling)
#        { 
#            sum( robustFct(((resp - multCurves(dose/scaleX, parm)) / scaleY)*weights), na.rm = rmNA)  
#        } else {
#        print(sum(robustFct(((resp / respScaling) - multCurves((dose / doseScaling), parm)) * weights), na.rm = rmNA))
#        print(c(parm, as.vector(multCurves((dose / doseScaling), parm)[1:5])))
#        print(as.vector(robustFct(((resp / respScaling) - multCurves((dose / doseScaling), parm)) * weights)))
#        print(multCurves((dose / doseScaling), parm))
        sum(robustFct(((resp / respScaling) - multCurves((dose / doseScaling), parm)) * weights), na.rm = rmNA)          
        # weights enter multiplicatively before being squared!
#        }
    }
    
    if (!is.null(dmf))
    {
        opdfct1 <- function(parm)
        {
#            apply(-2*(resp - multCurves(dose, parm))*dmf(dose, parm), 2, sum)
#            apply(-2*(resp - multCurves(dose, parm))*dmf(dose, parm), 2, appFct, cid)

            -2*((resp / respScaling) - multCurves((dose / doseScaling), parm)) * dmf((dose / doseScaling), parm)
        }
    } else {
        opdfct1 <- NULL
    }
       
    ## Defining additional self starter function (none needed)
    ssfct <- NULL

    ## Defining the log likelihood function
    llfct <- function(object)
    {
        degfre <- object$"sumList"$"lenData"  # "df.residual"  # object$summary[6]
        c( -(degfre/2)*(log(2*pi)+log(object$"fit"$"value")-log(degfre)+1), 
        object$"sumList"$"lenData" - object$"sumList"$"df.residual" + 1)
#        length(object$"fit"$"par") + 1 )
    }   
    
    ## Defining functions returning the residual variance, the variance-covariance matrix and the fixed effects estimates
    rvfct <- function(object)
    {
        object$"fit"$"value" / df.residual(object)  # object$"sumList"$"df.residual"
    }

    vcovfct <- function(object)
    {
#        scaledH <- (object$"fit"$"hessian")*(1/(2*object$"fit"$"ovalue"/object$"sumList"$"df.residual"))  # /2
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
#        solve((object$"fit"$"hessian")*(1/rvfct(object))/2)  
#        solve((object$"fit"$"hessian")*(1/(object$"fit"$"ovalue"/object$"sumList"$"df.residual"))/2)   
    }
    
    parmfct <- function(fit, fixed = TRUE)
    {
        fit$par
    }

    rstanfct <- function(object)
    {
        rep(1, object$"sumList"$"lenData")*sqrt(summary(object)$"resVar")
    }


    return(list(llfct = llfct, opfct = opfct, opdfct1 = opdfct1, ssfct = ssfct, 
    rvfct = rvfct, vcovfct = vcovfct, parmfct = parmfct, rstanfct = rstanfct))
}



