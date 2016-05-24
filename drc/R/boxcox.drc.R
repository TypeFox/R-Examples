"boxcox.drc" <- function(object, lambda = seq(-2, 2, by = 0.25), plotit = TRUE,
bcAdd = 0, method = c("ml", "anova"), level = 0.95, eps = 1/50, 
xlab = expression(lambda), ylab = "log-Likelihood", ...)
{
    method <- match.arg(method)
    
    ## Identifying the conditional or fixed-lambda approach
    if (identical(length(lambda), 1))
    {
        method <- "fixed"
    } 
    
    
    if (identical(method, "fixed"))
    {
        lv <- lambda
        ci <- c(NA, NA)
    }

    if (identical(method, "ml"))
    {
        ## Defining the likelihood function
        llFct <- function(object, lv) {
            yVec <- object$"data"[, 2]
            N <- length(yVec)  # df.residual(object)
            Ji <- yVec^(lv - 1)
            -N * log(sqrt(sum(residuals(object)^2)/N)) - N/2 + sum(log(Ji))
        }

        ## Fitting the model over a grid of lambda values
        lenlam <- length(lambda)
        llVec <- rep(NA, lenlam)
        for (i in 1:lenlam)
        {
            drcTemp <- try(update(object, bc = lambda[i], bcAdd = bcAdd), silent = TRUE)
            if (!inherits(drcTemp, "try-error")) 
            {
                llVec[i] <- llFct(drcTemp, lambda[i])  # logLik(drcTemp)
#                print(llVec[i])
            }
        }
        lv <- lambda[which.max(llVec)]
        ci <- boxcoxCI(lambda, llVec, level)    
#        llv <- max(llVec, na.rm = TRUE)

        ## Plotting the profile log-likelihood
        if (plotit)  # based on boxcox.default
        {
            plot(lambda, llVec, type ="l", xlab = xlab, ylab = ylab, ...)
        
            plims <- par("usr")
            y0 <- plims[3]
            llv <- max(llVec, na.rm = TRUE)
            lim <- llv - qchisq(level, 1)/2
        
            segments(lv, llv, lv, y0, lty=3)
            segments(ci[1], lim, ci[1], y0, lty = 3)  # lower limit
            segments(ci[2], lim, ci[2], y0, lty = 3)  # upper limit
        
            scal <- (1/10 * (plims[4] - y0))/par("pin")[2] 
            scx <- (1/10 * (plims[2] - plims[1]))/par("pin")[1] 
            text(lambda[1] + scx, lim + scal, " 95%") 
            abline(h = lim, lty = 3)
        } 
    }
     
    ## ANOVA-based approach 
    if (identical(method, "anova"))
    {
        dose <- object$dataList$"dose" 
        resp <- object$dataList$"resp"     
        curveid <- object$dataList$"curve"
        numCur <- length(unique(curveid)) 
     
        if (any(unlist(tapply(resp, dose, length)) < 2))
        {
            stop("ANOVA-based TBS approach requires replicates for each dose value")
        }
    
        ## Defining ANOVA formula
        afList <- anovaFormula(dose, resp, curveid, bcAdd)
        anovaForm <- afList$"anovaFormula"
        anovaData <- afList$"anovaData"

        profLik <- boxcox(anovaForm, lambda = lambda, plotit = plotit, data = anovaData)  
        # using boxcox in MASS         
                    
        lamVec <- profLik$x               
        llVec <- profLik$y
        lv <- lamVec[which.max(llVec)]
        ci <- boxcoxCI(lamVec, llVec, level)
    } 
     
    ## Updating the model fit
    retFit <- update(object, bc = lv, bcAdd = bcAdd)
    retFit$"boxcox" <- list(lambda = lv, ci = ci, bcAdd = bcAdd)
    retFit$call$bcVal <- lv
    retFit$call$bcAdd <- bcAdd    
#    retFit$boxcox[c(2, 3)] <- ci
    ## future: make boxcox and lambda into one component in the fit     
     
    ## Returning the result 
    invisible(retFit)
}

"boxcoxCI" <- 
function(x, y, level = 0.95)
{
    ## R lines taken from boxcox.default in the package MASS and then slightly modified
    xl <- x
    loglik <- y
    
    llnotna <- !is.na(loglik)
    xl <- xl[llnotna]
    loglik <- loglik[llnotna]
    
    m <- length(loglik)

    mx <- (1:m)[loglik == max(loglik)][1]
    Lmax <- loglik[mx]
    lim <- Lmax - qchisq(level, 1)/2

    ind <- range((1:m)[loglik > lim])

    xx <- rep(NA, 2)
    if(loglik[1] < lim) 
    {
        i <- ind[1]
        xx[1] <- xl[i - 1] + ((lim - loglik[i - 1]) *
                          (xl[i] - xl[i - 1]))/(loglik[i] - loglik[i - 1])

    }
    if(loglik[m] < lim) 
    {
        i <- ind[2] + 1
        xx[2] <- xl[i - 1] + ((lim - loglik[i - 1]) *
                          (xl[i] - xl[i - 1]))/(loglik[i] - loglik[i - 1])
    }
    return(xx)
}


## Defining ANOVA model formula
anovaFormula <- function(dose, resp, curveid, bcAdd)
{
    bcc <- rep(bcAdd, length(resp))    
    numCur <- length(unique(curveid))
    
    if (numCur > 1) 
    {
        anovaForm <- (resp + bcc) ~ offset(bcc) + factor(dose) * factor(curveid)
        alternative <- 2
    } else {
        anovaForm <- (resp + bcc) ~ offset(bcc) + factor(dose)
        alternative <- 1
    }
       
    list(anovaFormula = anovaForm, anovaData = data.frame(dose, resp, curveid, bcc))
}