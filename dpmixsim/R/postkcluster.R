##
## apply clustering to reduce number of total clusters 
##

postkcluster <-
function(mask, cx, clk=4, plot=TRUE)
{
    ncol <- ncol(mask)
    mask[mask != 0] <- cx
    cx <- as.vector(mask)
    ##-------------------------------------------------
    ## Option 1:  !!! For clustering with full matrix
    ## group=="clara"
    ## cat("\nOption 2: re-clustering using", clk, "medoids (cluster library)\n")
    clarax <- clara(cx, clk, metric="manhattan")  # clk:number of clusters
    cx <- clarax$clustering
    meds <- clarax$medoids
    tcx <- unique(cx)
    for(i in 1:length(tcx)) {
        nn <- which(cx == tcx[i])
        cx[nn] <- meds[i]
    }
    cat("colorcluster with",clk,"components:\n")
    print(table(cx))
    ## figname <- "test.eps"
    ##-------------------------------------------------
    ## Option 2:  !!! For clustering before mask
    ## tcx <- unique(cx)
    ##-------------------------------------------------
    if(plot) {
        cat("Display individual clusters: \n")
        tx <- table(cx)
        cat("palette=", as.integer(names(tx)),"\n")
        cl <- length(tx)     # !!! exclude background image if 0 ??
        ## cl <- length(tx) + 1
        op <-  options('warn')
        options('warn'=-1)
        cm <- matrix(1:cl, nrow=round(sqrt(cl)), byrow=TRUE)
        options('warn'=op$warn)
        ##---------------------
        par(ask=TRUE)
        par(mfcol=dim(cm), mar=c(0, 0, 0, 0) + 0.1)
        ## par(mfcol=dim(cm), mar=c(4, 3, 3, 2) + 0.1)
        ##---------------------
        ysim <- numeric(length(cx))
        for(i in 2:length(tx)) { # !!! exclude background image ??
        ## for(i in 1:length(tx)) 
            ysim[] <- 0
            vvi <- as.integer(names(tx[i]))
            nvi <- which(cx == vvi)
            ysim[nvi] <- cx[nvi]
            ims <- matrix(ysim, ncol=ncol)
            if(length(table(ims)) > 1){ # show only components with more than one color 
                image(ims, col=gray((0:255)/256), axes=FALSE)
            }
        }
        ## display regrouped classification
        ims <- matrix(cx, ncol=ncol)
        ## cat("display regrouped classification\n")
        rc <- range(ims)
        image(ims, col=gray((0:255)/256), zlim=c(0, rc[2]), axes=FALSE)
    }
    else {
        tx <- table(cx)
        cl <- length(tx) # !!! exclude background image if 0 ??
        ## cl <- length(tx) + 1
        op <-  options('warn')
        options('warn'=-1)
        cm <- matrix(1:cl, nrow=round(sqrt(cl)), byrow=TRUE)
        options('warn'=op$warn)
        ysim <- numeric(length(cx))
        ## for(i in 1:length(tx)) 
        for(i in 2:length(tx)) { # !!! exclude background image ??
            ysim[] <- 0
            vvi <- as.integer(names(tx[i]))
            nvi <- which(cx == vvi)
            ysim[nvi] <- cx[nvi]
            ims <- matrix(ysim, ncol=ncol)
        }
    }
    invisible(ims)
}

