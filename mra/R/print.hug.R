print.hug <- function( x, ... ){

nx <-   x$aux$nx 
ny <-   x$aux$ny

cap.coef <- round( x$capcoef, 5 )
se.cap <- round( x$se.capcoef, 5 )


if( ny > 0 ){
}

cat("Call:\n")
print(x$aux$call)
cat("\n")

if( ny <= 0 ){
    cat(" Capture and Recapture model:\n")
    coefmat <- paste( format( c(" Variable", names(cap.coef))), 
                      format( c(" Est", cap.coef) ),  
                      format( c(" SE", se.cap) ),
                      sep="  ")
} else {
    recap.coef <- cap.coef[ !x$remove ]
    se.recap <- se.cap[ !x$remove ]
    
    if( length(recap.coef) == 0 ){
        # All capture covars have been removed
        allgone <- TRUE
    } else {
        allgone <- FALSE
        nms <- paste("C:", names(recap.coef), sep="")
        fixed <- rep("(fixed)", length(recap.coef))
    }


    recap.coef <- c( recap.coef, round( x$recapcoef, 5 ))
    se.recap <- c( se.recap, round( x$se.recapcoef, 5 ))
    if( allgone ){
        nms <- names(x$recapcoef)
        fixed <- rep(NULL, length(recap.coef))
    } else {
        nms <- c( nms, paste("B:", names(x$recapcoef), sep=""))
        fixed <- c( fixed, rep("", length(x$recapcoef)))
    }    

    lc <- length(cap.coef)
    lr <- length(recap.coef)
    if( lc > lr ){
        recap.coef <- c(recap.coef, rep("", lc - lr))
        se.recap  <- c(se.recap, rep("", lc - lr))
        fixed <- c(fixed, rep("", lc - lr))
        nms <- c(nms, rep("", lc - lr))
    } else if (lc < lr) {
        cap.coef<- c(cap.coef, rep("", lr - lc))
        se.cap  <- c(se.cap, rep("", lr - lc))    
    }

    coefmat <- paste(   format( c(" Capture Model", names(cap.coef))), 
                        format( c(" Est", cap.coef) ),  
                        format( c(" SE", se.cap) ),
                        "  ",
                        format( c(" Recapture Model", nms)), 
                        format( c(" Est", recap.coef) ), 
                        format( c(" SE", se.recap) ),
                        sep="  ")
    coefmat <- paste(coefmat, c("",fixed), sep="")
                            
}
cat( paste(coefmat, "\n"))

cat("\nPopulation Size Estimate (se): ")
cat(paste( round(x$n.hat, 4), " (", round(x$se.n.hat,4), ")\n", sep=""))
cat(paste( x$n.hat.conf*100, "% confidence interval for population size: ", 
    round(x$n.hat.lower,2), " to ", round(x$n.hat.upper,2), "\n", sep=""))
cat(paste("Individuals observed: ", round(x$num.caught), "\n", sep=""))
cat(paste("Effective sample size: ", round(x$n.effective), "\n", sep=""))


cat(paste("\nMessage = ", x$message[2] ))
cat(paste("\nNumber of estimable coefficients (estimated) = ", x$df))
cat(paste("\nLog likelihood = ", x$loglike))
cat(paste("\nDeviance = ", x$dev))
cat(paste("\nAIC = ", x$aic))
cat(paste("\nAICc = ", x$aicc))
cat("\n")

invisible()

}
