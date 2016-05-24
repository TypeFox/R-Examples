plot.evmOpt <-
function(x, main=rep(NULL,4), xlab=rep(NULL,4), nsim=1000, alpha=.05, ...){
    if (!missing(main)){
        if (length(main) != 1 & length(main) != 4){
            stop("main should have length 1 or 4")
        }
        else if (length(main) == 1){ main <- rep(main, 4) }
    }

    if (all(sapply(x$data$D,ncol) == 1)){
        plot(ppevm(x, nsim=nsim, alpha=alpha), main=main[1], xlab=xlab[1])
        plot(qqevm(x, nsim=nsim, alpha=alpha), main=main[2], xlab=xlab[2])
        plotrl.evmOpt(x, main=main[3], xlab=xlab[3], smooth=FALSE, ...)
        plot(hist.evmOpt(x, main=main[4], xlab=xlab[4]))
    } else { # Covariates in the model
    
        np <- length(x$data$D)
        lp <- predict(x,type="lp", unique.=FALSE)$link
        Which <- as.logical(apply(lp[,1:np],2,var)) # identifies which cols have covariates

        x$data$y <- resid(x)
        x$threshold <- 0
        x$coefficients <- rep(0, length(x$data$D)) # phi not sigma, so 0 not 1
        plot(ppevm(x, nsim=nsim, alpha=alpha), main=main[1], xlab=xlab[1])
        plot(qqevm(x, nsim=nsim, alpha=alpha), main=main[2], xlab=xlab[2])

        
        for(i in (1:length(x$data$D))[Which]){
          ParName <- names(x$data$D[i])
          plot(lp[,i],resid(x),main=paste("Residuals vs fitted",ParName),xlab=paste("Fitted",ParName),ylab="Residuals")
          panel.smooth(lp[,i], resid(x), col.smooth=2)
        }
    }
    
    invisible()
}

