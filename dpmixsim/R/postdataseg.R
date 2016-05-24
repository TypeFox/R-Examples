##
## image segmentation
##
postdataseg <-
function(x, z, ngrid, dbg=FALSE)
{
    zsent <- c(z,-1) # add sentinel
    zjumps <- NULL
    cnt <- 1
    ## select sequential intervals with same zi
    for(k in 2:length(zsent)) {
        if(zsent[k] != zsent[k-1]) {
            zj <- c(zsent[k-1], cnt)
            zjumps <- rbind(zjumps, zj)
            cnt <- 1
        }
        else {
         cnt <- cnt+1
        }
    }
    row.names(zjumps) <- NULL    # to shutoff warning    
    zdf <- data.frame(zjumps, row.names=NULL)
    names(zdf) <- c("ix","freq")
    zc <- zdf$freq/ngrid
    zl <- c(0,cumsum(zc))
    if(dbg) {
        cat("thresholds:\n")
        print(zdf)
        print(zc)
        print(zl)
    }
    ##-----------------------------
    zt <- unique(z)
    ncolors <- length(zt)
    ## Option 0
    colors <- as.integer((1:(ncolors)) * 256/ncolors) # assigning colors ad-hoc
    ## colors <- as.integer((1:(ncolors)) * 256/(ncolors+1)) # assigning colors ad-hoc
    if(dbg) {
        cat("ncolors:", ncolors,"\n")
        cat("colors:", colors,"\n")
    }
    ##-----------------------------
    ## choose x values based on thresholds
    cx <- integer(length(x)) 
    for(k in 1:(length(zl)-1)) {
         j <- which(x >= zl[k] & x <= zl[k+1]) 
         cx[j] <- colors[zdf$ix[k]]
        if(dbg)
            cat("k:",k, "\tn.elements: ", length(j), "\tinterval: ", zl[k],":", zl[k+1],"\tcolor: ",
            colors[zdf$ix[k]],"\n")
    }
    ##
    ## im.classif <- matrix(cx,ncol=ncol)
    invisible(cx)
}

