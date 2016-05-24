
postimgclgrp <-
function(mask, cx, palcolor=TRUE)
{
    mask[mask != 0] <- cx
    cx <- mask
    ncol <- ncol(mask)
    cxm <- matrix(cx,ncol=ncol)
    rc <- range(cxm)
    par(ask=TRUE)
    if(!palcolor) {
        image(cxm, col=gray((0:rc[2])/rc[2]), zlim=c(0, rc[2]), axes=TRUE,
            main="clustered dpm classification")
    }
    else {
        ncolors <- length(table(cxm))
        image(cxm, col = 1:ncolors, zlim=c(0,256), axes=TRUE,
            main="clustered dpm classification")
    }

}

