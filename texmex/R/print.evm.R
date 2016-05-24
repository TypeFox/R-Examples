print.evmOpt <- function( x , digits=max(3, getOption("digits") - 3), ... ){
    cat( "Call: " )
    print( x$call, ... )
    print(x$family, verbose=FALSE)
    if ( is.null( x$penalty ) | x$penalty=="none" ){
        cat( "\nModel fit by maximum likelihood.\n" )
    }
    else {
        cat( "\nModel fit by penalized maximum likelihood.\n" )
    }
    if ( x$conv == 0 ) conv <- TRUE
    else conv <- FALSE
    cat( "\nConvergence:\t\t")
    cat(conv)

    if (x$rate < 1){
        cat( "\nThreshold:\t\t")
        cat(format(unname(x$threshold), digits=digits, ...))
        cat( "\nRate of excess:\t\t")
        cat(format(x$rate, digits=digits, ...))
    }

    cat("\n\n")
    
    if (x$penalty == "none"){
      wh <- t(format(c(x$loglik, AIC(x)), digits, ...))
      colnames(wh) <- c("Log. lik", "AIC")
      rownames(wh) <- ""
      print(wh, print.gap=2, quote=FALSE, justify="left")
    }
    else {
      wh <- t(format(c(x$loglik, x$ploglik, AIC(x)), digits, ...))
      colnames(wh) <- c("Log lik.", "Penalized log lik.", "AIC")
      rownames(wh) <- ""
      print(wh, print.gap=2, quote=FALSE, justify="left")
    }

    co <- cbind( coef(x), x$se )
    dimnames(co) <- list(names(coef(x)) , c("Value", "SE"))
    cat( "\n\nCoefficients:\n" )
    print.default(format(co, digits=digits, ...), print.gap=2, quote=FALSE)
    cat( "\n" )
    invisible()
}