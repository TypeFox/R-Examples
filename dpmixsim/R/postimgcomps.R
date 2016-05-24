
postimgcomps <-
function(mask, cx)
{
    cat("Display individual clusters: \n")
    mask[mask != 0] <- cx
    cx <- mask
    ncol <- ncol(mask)
    tx <- table(cx)
    cat("palette=", as.integer(names(tx)),"\n")
    cl <- length(tx)-1     # !!! exclude background image if 0 ??
    ## op <-  options('warn')
  	op <- par(no.readonly = TRUE) 
    options('warn'=-1)
    cm   <- matrix(1:cl, ncol=round(sqrt(cl)), byrow=TRUE)
    par(ask=TRUE)
    par(mfcol=dim(cm), mar=c(0, 0, 0, 0) + 0.1)
    ysim <- numeric(length(cx))
    for(i in 2:length(tx)) { # !!! exclude background image ??
        ysim[]    <- 0
        vvi       <- as.integer(names(tx[i]))
        nvi       <- which(cx == vvi)
        ysim[nvi] <- cx[nvi]
        ims       <- matrix(ysim, ncol=ncol)
        if(length(table(ims)) > 1){ # show only components with more than one color 
            image(ims, axes=FALSE, col=gray((0:2)/2))
        }
    }
    par(op)
}
