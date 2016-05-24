`print.ecespa.minconfit` <-
function (x, ...) {
    cat(paste("Minimum contrast fit ", "(", "object of class ", 
        dQuote("ecespa.minconfit"), ")", "\n", sep = ""))
    da <- x$dataname
    cat(paste("Model:PCP\n"))
    if(is.null(x$lambdaname)){
    cat(paste("Fitted by matching theoretical Kest function to", 
            x$dataname))
    }
    else cat(paste("Fitted by matching theoretical Kinhom function to", 
            x$dataname, "with lambda estimated by ", x$lambdaname))
    cat("\n")
    cat("Parameters fitted by minimum contrast ($sigma2 and $rho):\n")
    print(c(sigma2=x$sigma2, rho=x$rho,...))
    cat(paste("Domain of integration:", "[", signif(min(x$r), 
        4), ",", signif(max(x$r), 4), "]\n"))
    cat(paste("Exponents:", "p=", paste(signif(x$p, 3), ",", 
        sep = ""), "q=", signif(x$q, 3), "\n"))
    invisible(NULL)
}

