"drmEMeventtime" <- 
function(dose, resp, multCurves, doseScaling = 1)
{
    ## Defining the objective function                
    opfct <- function(c)  # dose, resp and weights are fixed
    {                      
        Fstart <- multCurves(dose[, 1] / doseScaling, c)
        dose2 <- dose[, 2]
#        IsFinite <- is.finite(dose2)
#        Fend <- rep(1, length(dose2))
#        Fend <- multCurves(dose[, 2] / doseScaling, c)
        Fend <- multCurves(dose2 / doseScaling, c)
        Fend[!is.finite(dose2)] <- 1 
        return( -sum(resp * log(Fend - Fstart)) )  
        # minus in front of sum() as maximization is done as minimization
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
        -object$"fit"$value,  # oops a constant is missing!
        object$"sumList"$"df.residual"
        )
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


"drmLOFeventtime" <- function()
{
    return(list(anovaTest = NULL, gofTest = NULL))
}
