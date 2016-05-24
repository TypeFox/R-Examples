`print.summary.splm` <-
function(x,digits= max(3, getOption("digits") - 2),width=getOption("width"),...) {


        cat(paste("Spatial panel",x$type,"model\n"))
        cat("\nCall:\n")
        print(x$call)
        cat("\nResiduals:\n")
        save.digits <- unlist(options(digits=digits))
        on.exit(options(digits=save.digits))
        print(sumres(x))

        if(!is.null(x$ErrCompTable)) {
            cat("\nError variance parameters:\n")
            printCoefmat(x$ErrCompTable,digits=digits,signif.legend=FALSE)
        }

        if(is.numeric(x$lambda)) {
            cat("\nEstimated spatial coefficient, variance components and theta:\n")
            print(x$lambda)
        }

        if(!is.null(x$ARCoefTable)) {
            cat("\nSpatial autoregressive coefficient:\n")
            printCoefmat(x$ARCoefTable,digits=digits,signif.legend=FALSE)
        }

        cat("\nCoefficients:\n")
        printCoefmat(x$CoefTable,digits=digits)
        cat("\n")

   

    invisible(x)
}

