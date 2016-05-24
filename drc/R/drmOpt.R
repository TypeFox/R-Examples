"drmOpt" <- 
function(opfct, opdfct1, startVec, optMethod, constrained, warnVal, 
upperLimits, lowerLimits, errorMessage, maxIt, relTol, opdfct2 = NULL, parmVec, traceVal, silentVal = TRUE,
matchCall)
## propagate "silentVal" from calling function? 
{
    ## Controlling the warnings
    options(warn = warnVal)   
    
    ## Calculating hessian
    if (is.null(opdfct2)) {hes <- TRUE} else {hes <- FALSE}

    ## Setting scaling parameters for optim()
    psVec <- abs(startVec)
    psVec[psVec < 1e-4] <- 1

    ## Derivatives are used
    {if (!is.null(opdfct1))
    {
        if (constrained)
        {
            nlsObj <- try(optim(startVec, opfct, opdfct1, hessian = hes, method = "L-BFGS-B", 
            lower = lowerLimits, upper = upperLimits, 
            control = list(maxit = maxIt, reltol = relTol, parscale = psVec)), silent = silentVal)
        } else {
            nlsObj <- try(optim(startVec, opfct, opdfct1, hessian = hes, method = optMethod, 
            control = list(maxit = maxIt, reltol = relTol, parscale = psVec)), silent = silentVal)
        }
        options(warn = 0)
        
        if (!inherits(nlsObj, "try-error")) 
        {
            nlsFit <- nlsObj
            nlsFit$convergence <- TRUE
        } else {
#            stop("Convergence failed")
            warning("Convergence failed. The model was not fitted!", call. = FALSE)

#            callDetail <- match.call()
#            if (is.null(callDetail$fct)) {callDetail$fct <- substitute(l4())}
            return(list(call = matchCall, parNames = parmVec, startVal = startVec, convergence = FALSE))
        }
        if (!hes) {nlsFit$hessian <- opdfct2(nlsFit$par)}

    ## Derivatives are not used
    } else {  

        if (constrained)
        {
#            print(lowerLimits)
#            print(upperLimits)
#            print(startVec)
#            print(opfct)
#            print(opfct(startVec))                    
            nlsObj <- try(optim(startVec, opfct, hessian = TRUE, method = "L-BFGS-B", 
            lower = lowerLimits, upper = upperLimits, 
            control = list(maxit = maxIt, parscale = psVec, reltol = relTol, trace = traceVal)), silent = silentVal)
            # parscale is needed for the example in methionine.Rd
        } else {
#            psVec <- abs(startVec)
#            psVec[psVec<1e-4] <- 1

            nlsObj <- try(optim(startVec, opfct, hessian = TRUE, method = optMethod, 
            control = list(maxit = maxIt, reltol = relTol, parscale = psVec, trace = traceVal)), silent = silentVal)

#            nlsObj0 <- try(optim(startVec, opfct, method=optMethod, 
#            control=list(maxit=maxIt, reltol=relTol, parscale=psVec)), silent=TRUE)
#            nlsObj <- try(optim(nlsObj0$par, opfct, hessian=TRUE, method=optMethod, 
#            control=list(maxit=maxIt, reltol=relTol)), silent=TRUE)
        }
        options(warn = 0)
        
        if (!inherits(nlsObj, "try-error")) 
        {
            nlsFit <- nlsObj
            nlsFit$convergence <- TRUE            
        } else {  # to avoid an error if used in a loop
            if (errorMessage) 
            {
                stop("Convergence failed")
            } else {
                warning("Convergence failed. The model was not fitted!", call. = FALSE)
            }

#            callDetail <- match.call()
#            if (is.null(callDetail$fct)) {callDetail$fct <- substitute(LL.4())}
            return(list(call = matchCall, parNames = parmVec, startVal = startVec, convergence = FALSE))
        }
    }}
    
#    nlsFit$ofvalue <- nlsFit$value
    nlsFit$ovalue <- nlsFit$value  # used in the var-cov matrix ... check
#    nlsFit$value <- opfct(nlsFit$par, scaling = FALSE)  # used in the residual variance ... check    
    nlsFit$value <- opfct(nlsFit$par)

    ## Returning the fit
    return(nlsFit)
}
