## frair methods

# Print for objects of class frfit
print.frfit <- function(x, ...){
    cat('FUNCTIONAL RESPONSE FIT\n')
    cat(paste0('\nResponse:            ', x$response))
    cat(paste0('\nDescription:         ', as.character(frair_responses(show=FALSE)[[x$response]][2])))
    cat(paste0('\nOptimised variables: ', paste0(x$optimvars, collapse=', ')))
    cat(paste0('\nFixed variables:     ', ifelse(test=!is.null(x$fixedvars), yes=paste(x$fixedvars, collapse=', '), no='NA')))
    cat('\n')
    cat('\nCoefficients:\n')
    print(round(x$coefficients, 3))
    cat('\nNOTE: It is recomended you inspect the raw fit too (see: ?frair_fit)\n')
}

plot.frfit <- function(x, xlab=x$xvar, ylab=x$yvar, ...){
    plot(x$x, x$y, xlab=xlab, ylab=ylab, ...)
}

lines.frfit <- function(x, tozero=FALSE, ...){
    fitfun <- get(x$response, pos = "package:frair")
    if(tozero){
        zero_answer <- fitfun(0, as.list(x$coefficients))
        if(is.na(zero_answer)){
            warning(c("The supplied function is undefined at zero.\n",
                      "   Plotting to a minimum of 1e-04 instead."))
            lowval <- 1e-04
        } else {
            lowval <- 0
        }
        newx <- seq(from=lowval, to=max(x$x), length.out = 50)
    } else {
        newx <- seq(from=min(x$x), to=max(x$x), length.out = 50)
    }
    newy <- fitfun(newx, as.list(x$coefficients))
    lines(x=newx, y=newy, ...)
}