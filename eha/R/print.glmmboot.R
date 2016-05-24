print.glmmboot <- function(x,
                         digits = max(3, getOption("digits") - 3),
                         na.print = "",
                         ...){ 
    
    cat("\nCall: ", deparse(x$call), "\n\n")
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    if (length(x$coefficients)){
        coef <- x$coefficients
        se <- x$sd
        tmp <- cbind(coef,
                     se,
                     coef/se,
                     signif(1 - pchisq((coef/se)^2, 1), digits - 1)
                     )
        dimnames(tmp) <- list(names(coef),
                              c("coef", "se(coef)", "z", "Pr(>|z|)")
                              )
        cat("\n")
        prmatrix(tmp)
    }
    
    if (x$boot_rep){
        cat("\n Bootstrap p-value for fixed mixing: ",
            x$bootP, "(", x$boot_rep, ")\n")
    }
    
    cat("\nResidual deviance:",
        format(signif(x$deviance, digits)), "on",
        x$df.residual, "degrees of freedom", 
        "\tAIC:",
        format(signif(x$aic, digits)), "\n")
}
