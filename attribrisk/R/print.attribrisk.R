print.attribrisk <- function(x, ...) {
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
    }

    if (inherits(x, "coxph"))
        cat("\nn = ", length(x$fit$residuals), "\n")
    else if (inherits(x, 'glm'))
        cat('\n Degrees of freedom:', x$fit$df.null, '\n')

    # Put together a print vector: a 1 row matrix with names
    if (is.null(x$var)) { # rare case -- no variance done
        cat("Attributable risk:", format(x$attribrisk, ...), "\n")
    }
    else {
        temp <- c(x$attribrisk, sqrt(x$var))
        if (!is.null(x$boot.ci)) {
            # This next part is dangerous in a way.  We want to print
            #  the boot.ci results, but without the "call" part since
            #  our internal call is irrelevant.  If the format of boot.ci
            #  objects were to change the below would break.
            #
            # grab the first type found, if there are multiple
            type <- match(names(x$boot.ci), c("normal", "basic", "student",
                                            "percent", "bca"), nomatch=0)
            type <- min(which(type>0))
            temp <- c(temp, x$boot.ci[[type]][4:5])
        }
        else {
            z <- qnorm((1-x$conf)/2)
            temp <- c(temp, temp[1] + c(-1,1)*z*temp[2])
        }

        temp <- matrix(temp, nrow=1)
        
        dimnames(temp) <- list("attributable risk", 
                               c("coefficient", "std. err", 
                                 paste(c("lower", "upper"), format(x$conf))))
        print(temp, ...)
    }
    invisible(x)
}


