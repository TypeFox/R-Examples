"drmEMPoisson" <- 
function(dose, resp, multCurves, startVec, weightsVec, doseScaling = 1)
{

    ## Finding indices for doses that give contribution to likelihood function
#    iv <- ( (multCurves(dose, startVec) > zeroTol) & (multCurves(dose, startVec) < 1-zeroTol) )


    ## Defining the objective function                
    opfct <- function(c)  # dose, resp and weights are fixed
    {                      
        lambda <- weightsVec * multCurves(dose / doseScaling, c)
        return( -sum(-lambda + resp*log(lambda)))
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
    
       
    ## Defining functions returning the residual variance, the variance-covariance and the fixed effects estimates
    rvfct <- NULL

    vcovfct <- function(object)
    {
        solve(object$fit$hessian)    
    }
    
    parmfct <- function(fit, fixed = TRUE)
    {
        fit$par
    }


    ## Returning list of functions
    return(list(llfct = llfct, opfct = opfct, ssfct = ssfct, rvfct = rvfct, vcovfct = vcovfct, 
    parmfct = parmfct))
}


"drmLOFPoisson" <- function()
{
    return(list(anovaTest = NULL, gofTest = NULL))
}
