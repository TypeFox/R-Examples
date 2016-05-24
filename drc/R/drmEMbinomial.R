"drmEMbinomial" <- 
function(dose, resp, multCurves, startVec, robustFct, weights, rmNA, zeroTol = 1e-12, 
doseScaling = 1, respScaling = 1)
{
    ## Finding indices for doses that give contribution to likelihood function
    iv <- ( (multCurves(dose/doseScaling, startVec) > zeroTol) & (multCurves(dose/doseScaling, startVec) < 1-zeroTol) )

    ## Defining the objective function                
    opfct <- function(c)  # dose, resp and weights are fixed
    {                      
#        prob <- (multCurves(dose / doseScaling, c))[iv]
#        prob <- multCurves(dose, c)

#         prob0 <- (multCurves(dose / doseScaling, c))        
#         iv <- (prob0 > zeroTol) & (prob0 < (1 - zeroTol))
#         prob <- prob0[iv]

#         print(log(prob/(1-prob[]))) 
#        print(-sum((resp*weights)[iv]*log(prob/(1-prob))+weights[iv]*log(1-prob)))
#        return( -sum((resp2*weights2)*log(prob/(1-prob))+weights2*log(1-prob)) )
        
        prob <- multCurves(dose / doseScaling, c)
        omZT <- 1 - zeroTol
        prob[prob > omZT] <- omZT
        prob[prob < zeroTol] <- zeroTol
        -sum((resp * weights) * log(prob / (1 - prob)) + (weights * log(1 - prob)))        
#        -sum((resp * weights)[iv] * log(prob / (1 - prob)) + (weights[iv] * log(1 - prob)))
        
#        print(sum(retVec[is.finite(retVec)])) 
#        sum(retVec[is.finite(retVec)])
#        print(retVec)
#        sum(retVec)

#        return( -sum((resp*weights)[iv]*log(prob/(1-prob))+weights[iv]*log(1-prob)) )
    }    
    
    ## Defining self starter function
    ssfct <- NULL

    ## Defining the log likelihood function
    llfct <- function(object)
    {
#        total <- (object$"data")[iv, 5]
#        success <- total*(object$"data")[iv, 2]    
        total <- (object$"data")[, 5]
        success <- total*(object$"data")[, 2]    

        
        c(sum(log(choose(total, success))) - object$"fit"$"ovalue",  # object$"fit"$"ofvalue", 
        object$"sumList"$"lenData" - df.residual(object))  # object$"sumList"$"df.residual")
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

#
#    ## Modifying ANOVA test (removing dose=0 and dose=Inf)
#    anovaTest2 <- function(formula, ds) {anovaTest(formula, ds[iv, ])}


    ## Returning list of functions
    return(list(llfct = llfct, opfct = opfct, ssfct = ssfct, rvfct = rvfct, 
    vcovfct = vcovfct, parmfct = parmfct))  # , anovaTest2=anovaTest2))
}


