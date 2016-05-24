plotinteractpreds <- function(out, seplot=TRUE, ylim=NULL, main=NULL, xlab=NULL, ylab=NULL, legend=TRUE, placement="bottomright", lwd=3, add=FALSE, addat=FALSE, mfrow=NULL, linecol=NULL, secol=NULL, showbynamelegend=FALSE, showatnamelegend=FALSE, ...){
    oldmfrow <- par()$mfrow
    if(class(out)!="interactpreds")
        warning("This function may not work with data that is not generated using the findwtdinteraction function")
    if(!is.null(mfrow))
        par(mfrow=mfrow)
    dvname <- out$Meta$dvname
    across <- out$Meta$across
    by <- out$Meta$by
    at <- out$Meta$at
    mins <- lapply(1:length(out$Means), function(x) out$Means[[x]]-1.96*out$SEs[[x]])
    maxs <- lapply(1:length(out$Means), function(x) out$Means[[x]]+1.96*out$SEs[[x]])
    if(is.null(ylim))
        ylim <- c(min(unlist(mins), na.rm=TRUE), max(unlist(maxs), na.rm=TRUE))
    atlevs <- try(names(out$Means))
    bylevs <- try(rownames(out$Means[[1]]))
    acrosslevs <- try(names(out$Means[[1]]))
    if(!is.null(bylevs)){
        acrosslevs <- try(colnames(out$Means[[1]]))
        hasby=TRUE
        bylevnames <- bylevs
    }
    if(is.null(bylevs)){
        hasby <- FALSE
        bylevs <- "All"
    }
    accnumeric <- suppressWarnings(sum(as.character(as.numeric(acrosslevs))==acrosslevs, na.rm=TRUE)==length(acrosslevs))
    if(accnumeric==TRUE){
        acrosslevs <- as.numeric(acrosslevs)
        acrossvals <- as.numeric(acrosslevs)
    }
    if(accnumeric==FALSE)
        acrossvals <- 1:length(acrosslevs)
    if(is.null(xlim))
        xlim <- c(min(acrossvals, na.rm=TRUE), max(acrossvals, na.rm=TRUE))
    hasat <- !(is.null(atlevs) || length(atlevs)==1)
    if(is.null(linecol) & addat==TRUE)
        linecol <- gray((1:length(atlevs))/(2*length(atlevs)))
    if(is.null(linecol))
        linecol <- "black"
    if(is.null(secol) & addat==TRUE)
        secol <- gray(.25+(1:length(atlevs))/(2*length(atlevs)))
    atlegend <- atlevs
    bylegend <- bylevs
    if(showatnamelegend==TRUE)
        atlegend <- paste(atlevs, at)
    if(showbynamelegend==TRUE)
        bylegend <- paste(bylevs, by)
    premain <- main
    for(a in 1:length(atlevs)){
        if(is.null(premain)){
            mainp <- paste("Interaction Plot of", dvname, "Across Levels of\n", across)
            if(!is.null(bylevs) & length(bylevs)>1)
                mainp <- paste(mainp, "By", by)
            if(hasat==TRUE & addat==FALSE)
                main <- paste(mainp, "At", at, "=", atlevs[a])
            if(hasat==TRUE & addat==TRUE)
                main <- paste(mainp, "At Each Level Of", at)
            if(hasat==FALSE)
                main <- mainp
        }
        if(is.null(ylab))
            ylab <- dvname
        if(is.null(xlab))
            xlab <- across
        if(add==FALSE){
            plot(acrossvals, acrossvals, type="n", ylim=ylim, main=main, ylab=ylab, xlab=xlab, axes=FALSE)#, ...)
            axis(1, at=c(-999,999))
            axis(1, at=acrossvals, labels=acrosslevs)
            axis(2, at=c(-999,999))
            axis(2)
            axis(3, at=c(-999,999))
            axis(4, at=c(-999,999))
        }
        if(seplot==TRUE){
            if(hasby==TRUE){
                for(i in 1:length(bylevs)){
                    if(seplot==TRUE & hasat==FALSE)
                        polygon(c(acrossvals, rev(acrossvals)), c(maxs[[a]][i,], rev(mins[[a]][i,])), density=30, angle=45+10*i, col=secol)
                    if(seplot==TRUE & hasat==TRUE & addat==FALSE)
                        polygon(c(acrossvals, rev(acrossvals)), c(maxs[[a]][i,], rev(mins[[a]][i,])), density=30, angle=45+10*i, col=secol)
                    if(seplot==TRUE & hasat==TRUE & addat==TRUE)
                        polygon(c(acrossvals, rev(acrossvals)), c(maxs[[a]][i,], rev(mins[[a]][i,])), density=30, angle=45*a+10*i, col=secol[a])
                }
            }
            if(hasby==FALSE)
                polygon(c(acrossvals, rev(acrossvals)), c(maxs[[a]], rev(mins[[a]])), density=30, angle=45*a, col=secol[a])
        }
        if(addat==TRUE)
            add <- TRUE
        if(length(linecol)==1 & hasat==TRUE & addat==TRUE)
            linecol <- rep(linecol, length(atlevs))
        if(length(linecol)==1)   
            linecol <- rep(linecol, length(bylevs))
        if(hasby==TRUE){
            for(i in 1:length(bylevs)){
                if(hasat==FALSE || (hasat==TRUE & addat==FALSE))
                    lines(acrossvals, out$Means[[a]][i,], lty=i, lwd=lwd, col=linecol[i])
                if(hasat==TRUE & addat==TRUE)
                    lines(acrossvals, out$Means[[a]][i,], lty=i, lwd=lwd, col=linecol[a])            
                if(legend==TRUE){
                    if(hasat==FALSE)
                        legend(x=placement, legend=bylegend, lty=1:length(bylevs), lwd=lwd, col=linecol)
                    if(hasat==TRUE & addat==FALSE)
                        legend(x=placement, legend=bylegend, lty=1:length(bylevs), lwd=lwd, col=linecol)
                }
            }
        }
        if(hasby==FALSE){
            if(hasat==FALSE || (hasat==TRUE & addat==FALSE))
                lines(acrossvals, out$Means[[a]], lty=1, lwd=lwd, col=linecol)
            if(hasat==TRUE & addat==TRUE)
                lines(acrossvals, out$Means[[a]], lty=a, lwd=lwd, col=linecol[a])
        }
    }
    if(legend==TRUE & hasat==TRUE & addat==TRUE & hasby==TRUE)
        legend(x=placement, legend=c(bylegend, atlegend), lty=c(1:length(bylevs), rep(1, length(atlevs))), lwd=lwd, col=c(rep("black", length(bylevnames)), linecol))
    if(legend==TRUE & hasat==TRUE & addat==TRUE & hasby==FALSE)
        legend(x=placement, legend=c(atlegend), lty=1:length(atlevs), lwd=lwd, col=linecol)
    if(!is.null(mfrow))
        par(mfrow=oldmfrow)
}
