gam.style <-
function(x, parms, column, baseline=mean(x[,column]),
         epsilon=1e-5, seg.len=0.02, seg.cols="black", plot=TRUE,
         return.results=FALSE, trim=0, ...)
{
    effect <-
    function(x, column, parms, baseline, ensemble.average, trim)
    {
        parms$lower <- -Inf
        x.baseline <- x
        x.baseline[,column] <- baseline
        p.baseline <- qrnn.predict(x.baseline, parms)
        p <- qrnn.predict(x, parms)
        if(ensemble.average){
            p.baseline <- apply(p.baseline, 1, censored.mean,
                                lower=parms$lower, trim=trim)
            p <- apply(p, 1, censored.mean, lower=parms$lower, trim=trim)
        }
        p - p.baseline
    }
    partial <-
    function(x, column, parms, epsilon, ensemble.average, trim)
    {
        parms$lower <- -Inf
        x.plus <- x.minus <- x
        x.plus[,column] <- x.plus[,column] + epsilon
        x.minus[,column] <- x.minus[,column] - epsilon
        p.plus <- qrnn.predict(x.plus, parms)
        p.minus <- qrnn.predict(x.minus, parms)
        if(ensemble.average){
            p.plus <- apply(p.plus, 1, censored.mean, lower=parms$lower,
                            trim=trim)
            p.minus <- apply(p.minus, 1, censored.mean, lower=parms$lower,
                             trim=trim)
        }
        (p.plus - p.minus)/(2*epsilon)
    }
    effects <- effect(x, column, parms, baseline, plot, trim)
    partials <- partial(x, column, parms, epsilon, plot, trim)
    if(plot){
        if(is.null(colnames(x)))
            colnames(x) <- paste("x", 1:ncol(x), sep="")
        xlab <- colnames(x)[column]
        x.var <- x[,column]
        if(length(seg.cols)==1) seg.cols <- rep(seg.cols, length(x.var))
        ylab <- paste("Effects: tau =", parms$tau)
        theta <- atan(partials)
        ymin <- min(effects)
        ymax <- max(effects)
        xmin <- min(x.var)
        xmax <- max(x.var)
        aspect <- (ymax - ymin)/(xmax - xmin)
        xdev <- seg.len*(xmax-xmin)*cos(theta)
        ydev <- seg.len*(xmax-xmin)*sin(theta)
        scale <- sqrt(xdev**2 + (ydev/aspect)**2)
        xdev <- xdev*(seg.len*(xmax-xmin))/scale
        ydev <- ydev*(seg.len*(xmax-xmin))/scale
        dev.new()
        plot(x.var, effects, type="n", xlab=xlab, ylab=ylab, ...)
        for(case in seq_along(x.var)){
            xi <- x.var[case]
            yi <- effects[case]
            xd <- xdev[case]
            yd <- ydev[case]
            segments(xi, yi, xi+xd, yi+yd, col=seg.cols[case])
        }
        abline(v=baseline, lty=3)
    }
    if(return.results) return(list(effects=effects, partials=partials))
}
