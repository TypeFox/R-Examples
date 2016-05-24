gam.style <-
function(x, weights, column, baseline=mean(x[,column]),
         epsilon=1e-5, seg.len=0.02, seg.cols="black", plot=TRUE,
         return.results=FALSE, ...)
{
    effect <-
    function(x, column, weights, baseline)
    {
        x.baseline <- x
        x.baseline[,column] <- baseline
        p.baseline <- monmlp.predict(x.baseline, weights)
        p <- monmlp.predict(x, weights)
        effs <- p - p.baseline
        if(!is.null(attr(p.baseline, "ensemble"))){
            p.ens.baseline <- attr(p.baseline, "ensemble")
            p.ens <- attr(p, "ensemble")
            effs.ens <- list()
            for(i in 1:length(p.ens.baseline)){
                effs.ens[[i]] <- p.ens[[i]] - p.ens.baseline[[i]]
            }
        attr(effs, "ensemble") <- effs.ens
        }
        effs
    }
    partial <-
    function(x, column, weights, epsilon)
    {
        x.plus <- x.minus <- x
        x.plus[,column] <- x.plus[,column] + epsilon
        x.minus[,column] <- x.minus[,column] - epsilon
        p.plus <- monmlp.predict(x.plus, weights)
        p.minus <- monmlp.predict(x.minus, weights)
        pds <- (p.plus - p.minus)/(2*epsilon)
        if(!is.null(attr(p.plus, "ensemble"))){
            p.ens.plus <- attr(p.plus, "ensemble")
            p.ens.minus <- attr(p.minus, "ensemble")
            pds.ens <- list()
            for(i in 1:length(p.ens.plus)){
                pds.ens[[i]] <- (p.ens.plus[[i]] - p.ens.minus[[i]])/(2*epsilon)
            }
        attr(pds, "ensemble") <- pds.ens
        }
        pds
    }
    effects <- effect(x, column, weights, baseline)
    partials <- partial(x, column, weights, epsilon)
    if(plot){
        x.var <- x[,column]
        if(is.null(colnames(x)))
            colnames(x) <- paste("x", 1:ncol(x), sep="")
        xlab <- colnames(x)[column]
        if(length(seg.cols)==1) seg.cols <- rep(seg.cols, length(x.var))
        for(predictand in 1:ncol(effects)){
            ylab <- paste("Effects: predictand", predictand)
            theta <- atan(partials[,predictand])
            ymin <- min(effects[,predictand])
            ymax <- max(effects[,predictand])
            xmin <- min(x.var)
            xmax <- max(x.var)
            aspect <- (ymax - ymin)/(xmax - xmin)
            xdev <- seg.len*(xmax-xmin)*cos(theta)
            ydev <- seg.len*(xmax-xmin)*sin(theta)
            scale <- sqrt(xdev**2 + (ydev/aspect)**2)
            xdev <- xdev*(seg.len*(xmax-xmin))/scale
            ydev <- ydev*(seg.len*(xmax-xmin))/scale
            dev.new()
            plot(x.var, effects[,predictand], type="n", xlab=xlab, ylab=ylab,
                 ...)
            for(case in seq_along(x.var)){
                xi <- x.var[case]
                yi <- effects[case,predictand]
                xd <- xdev[case]
                yd <- ydev[case]
                segments(xi, yi, xi+xd, yi+yd, col=seg.cols[case])
            }
            abline(v=baseline, lty=3)
        }
    }
    if(return.results) return(list(effects=effects, partials=partials))
}
