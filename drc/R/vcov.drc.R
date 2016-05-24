"vcov.drc" <-
function(object, ..., corr = FALSE, od = FALSE, pool = TRUE, unscaled = FALSE)
{
    ## Defining function for calculating variance-covariance matrix
    contData <- identical(object$"type", "continuous") || identical(object$"type", "standard")
    if (contData)
    {
        vcovfct <- vcCont
    } else {
        vcovfct <- vcDisc
    }

    ## Retrieving the estimated variance-covariance matrix for the parameter estimates
    if (!corr)
    {
#        summary(object)$"varMat"
        if (!is.null(object$"objList"))
        {
#            require(magic, quietly = TRUE)

            if ((contData) && (pool))
            {
                vcovfct <- function(x){vcCont(x) / (2 * rse(x, TRUE))}  
                # no individual scaling, only for continuous data

                vcMat <- do.call("adiag", lapply(object$"objList", vcovfct))
                if (!unscaled)
                {
                    vcMat <- vcMat * (2 * (object$"minval" / df.residual(object)))
                    # scaling based on all fits
                }
            } else {
                vcMat <- do.call("adiag", lapply(object$"objList", vcovfct))
            }
#            do.call("adiag", lapply(object$"objList", object$"estMethod"$"vcovfct"))
#            vcMat <- do.call("adiag", lapply(object$"objList", vcovfct))
#            if (contPool)
#            {
#                vcMat <- vcMat * (2 * (object$"minval" / df.residual(object)))
#                # scaling based on all fits
#            }        
            return(vcMat)     
        } else {
            if ((contData) && (unscaled))
            {
                return(vcovfct(object) / (2 * rse(object, TRUE)))
            } else {
                varMat <- vcovfct(object)
                ## Adjusting for over-dispersion using the Pearson statistic
                if ((identical(object$"type", "binomial")) && (od)) 
                { 
                    gof <- drmLOFbinomial()$"gofTest"(object$"dataList"$"resp", 
                    weights(object), fitted(object), df.residual(object))
                
                    varMat <- varMat * (gof[1] / gof[2])
                }
                return(varMat)
            }
        }
    } else { ## Calculating correlation matrix
        corrFct <- function(object)
        {
            vcMat <- (object$"estMethod")$"vcovfct"(object)
            diage <- sqrt(diag(vcMat))
            vcMat / (outer(diage, diage))
        }
        if (!is.null(object$"objList"))
        {
#            require(magic, quietly = TRUE)        
            do.call("adiag", lapply(object$"objList", corrFct)) 
        } else {
            corrFct(object)
        }       
    }    
}

"vcCont" <- function(object)
{
#    scaledH <- (object$"fit"$"hessian") / (2 * rvfct(object))
    scaledH <- (object$"fit"$"hessian") / (2 * rse(object, TRUE))
    invMat <- try(solve(scaledH), silent = TRUE)
    
    if (inherits(invMat, "try-error"))
    {
#        cat("Note: Variance-covariance matrix regularized\n")
        ## More stable than 'solve' (suggested by Nicholas Lewin-Koh - 2007-02-12)
        ch <- try(chol(scaledH), silent = TRUE)  
        ## "silent" argument added after report by Xuesong Yu - 2010-03-09
        if (inherits(ch, "try-error")) 
        {
#            ch <- try(chol(0.99 * object$fit$hessian + 0.01 * diag(dim(object$fit$hessian)[1])), silent = TRUE)
            ch <- try(chol(0.99 * scaledH + 0.01 * diag(dim(scaledH)[1])), silent = TRUE)  # 2012-06-22
        }
        ## Try regularizing if the varcov is unstable
        if (!inherits(ch, "try-error")) 
        {
            return(chol2inv(ch))
        } else {
            numRows <- dim(scaledH)[1]
            return(matrix(NA, numRows, numRows))
        }
    } else {
        return(invMat)
    }
}

"vcDisc" <- function(object)
{
    solve(object$fit$hessian)    
}