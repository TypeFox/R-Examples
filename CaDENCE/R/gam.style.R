gam.style <-
function(x, fit, column, baseline=mean(x[,column]), additive.scale=FALSE,
         epsilon=1e-5, seg.len=0.02, seg.cols="black", plot=TRUE,
         return.results=FALSE, ...)
{
    effect <-
    function(x, column, fit, baseline)
    {
        x.baseline <- x
        x.baseline[,column] <- baseline
        p.baseline <- cadence.predict(x.baseline, fit)
        p <- cadence.predict(x, fit)
        p - p.baseline
    }
    partial <-
    function(x, column, fit, epsilon)
    {
        x.plus <- x.minus <- x
        x.plus[,column] <- x.plus[,column] + epsilon
        x.minus[,column] <- x.minus[,column] - epsilon
        p.plus <- cadence.predict(x.plus, fit)
        p.minus <- cadence.predict(x.minus, fit)
        (p.plus - p.minus)/(2*epsilon)
    }
    if (!("W1" %in% names(fit)))
        stop("\"fit\" must be a single list element from \"cadence.fit\"")
    parameters <- attr(fit, "distribution")$parameters
    if (additive.scale) {
        for(param in 1:length(parameters)) {
            attr(fit, "distribution")$output.fcns[[param]] <- identity
        }
    }   
    effects <- effect(x, column, fit, baseline)
    partials <- partial(x, column, fit, epsilon)
    colnames(effects) <- colnames(partials) <- parameters
    if(plot){
        if(is.null(colnames(x)))
            colnames(x) <- paste("x", 1:ncol(x), sep="")
        xlab <- colnames(x)[column]
        x.var <- x[,column]
        if(length(seg.cols)==1) seg.cols <- rep(seg.cols, length(x.var))
        for(param in 1:ncol(effects)){
            ylab <- paste("Effects:", parameters[param])
            theta <- atan(partials[,param])
            ymin <- min(effects[,param])
            ymax <- max(effects[,param])
            xmin <- min(x.var)
            xmax <- max(x.var)
            aspect <- (ymax - ymin)/(xmax - xmin)
            xdev <- seg.len*(xmax-xmin)*cos(theta)
            ydev <- seg.len*(xmax-xmin)*sin(theta)
            scale <- sqrt(xdev**2 + (ydev/aspect)**2)
            xdev <- xdev*(seg.len*(xmax-xmin))/scale
            ydev <- ydev*(seg.len*(xmax-xmin))/scale
            dev.new()
            plot(x.var, effects[,param], type="n", xlab=xlab, ylab=ylab, ...)
            for(case in seq_along(x.var)){
                xi <- x.var[case]
                yi <- effects[case,param]
                xd <- xdev[case]
                yd <- ydev[case]
                segments(xi, yi, xi+xd, yi+yd, col=seg.cols[case])
            }
            abline(v=baseline, lty=3)
        }
    }
    if(return.results) return(list(effects=effects, partials=partials))
}
