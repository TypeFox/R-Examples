###############################################################################
## Computation of distances for RFLP data
###############################################################################

## x: data.frame with RFLP data
## distfun: function to compute distance (cf. ?dist)
RFLPplot <- function(x, nrBands, nrMissing, distfun = dist, 
		     hclust.method = "complete", 
                     mar.bottom = 5, cex.axis = 0.5,
		     colBands, xlab = "", ylab = "molecular weight", 
		     ylim, ...){
    stopifnot(is.data.frame(x))
    stopifnot(is.function(distfun))
    
    if(missing(nrBands))
        stop("Number of Bands 'nrBands' is missing.")
        
    x1 <- split(x, x$Sample)
    x1.bands <- sapply(x1, nrow)

    if(missing(nrMissing)){
        ind <- x1.bands == nrBands
        dend.ord <- order.dendrogram(as.dendrogram(hclust(RFLPdist(x, distfun = distfun, 
                                                                   nrBands = nrBands), 
                                                          method = hclust.method)))
    }else{
        ind <- x1.bands %in% c(nrBands:(nrBands+nrMissing))
        dend.ord <- order.dendrogram(as.dendrogram(hclust(RFLPdist2(x, distfun = distfun, 
                                                                   nrBands = nrBands,
                                                                   nrMissing = nrMissing), 
                                                          method = hclust.method)))
    }
    temp <- x1[ind]
    temp1 <- do.call("rbind", temp)
    
    par(mar = c(mar.bottom, 4, 4, 2) + 0.1)    
    if(missing(ylim)){
        lo <- 10*trunc(min(temp1$MW)/10)
        up <- 10*ceiling(max(temp1$MW)/10)
        if(up-lo < 1) up <- up + 1
        ylim <- c(lo, up)
    }
    
    plot(NA, NA, xlim = c(0.25,length(temp)+0.75), ylim = ylim,
         xlab = xlab, ylab = ylab, axes = FALSE, ...)
    At <- round(seq(ylim[1], ylim[2], length = 11), 0)
    axis(2, at = At, labels = as.character(At), las = 2)
    axis(1, at = 1:length(temp), labels = names(temp)[dend.ord], las = 2, 
         cex.axis = cex.axis)
    abline(h = At, col = "grey", lty = 2)
    box()
    
    if(missing(nrMissing)){
        title(paste("Samples with", nrBands, "bands"))
    }else{
        if(nrMissing == 1){
            title(paste("Samples with", nrBands, "or", nrBands+1, "bands"))
        }else{
            title(paste("Samples with", nrBands, ", ..., ", nrBands+nrMissing, "bands"))
        }
    }
    rm(temp1)
    if(missing(colBands)){
	if(length(temp) <= 9){
	    mycol <- brewer.pal(max(3, length(temp)), "Set1")
	}else{
	    mycol1 <- brewer.pal(9, "Set1")
	    mycol <- colorRampPalette(mycol1)(length(temp))
	}
    }else{
	if((length(colBands) != 1) && (length(colBands) != length(temp)))
	    stop("Length of 'colBands' is", length(colBands), "and should be 1 or", length(temp))
	if(length(colBands) == 1){ 
	    mycol <- rep(colBands, length(temp))
        }else{
            mycol <- colBands
        }
    }

    temp.o <- temp[dend.ord]
    for(i in seq_along(temp.o)){
        temp.mw <- temp.o[[i]]$MW
        reps <- length(temp.mw)
        matlines(rbind(rep(i-0.25, reps), rep(i+0.25, reps)), 
                 rbind(temp.mw, temp.mw), 
                 lwd = 3, lty = 1, col = mycol[i])
    }
    invisible()
}
