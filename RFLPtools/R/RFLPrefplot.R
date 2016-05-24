###############################################################################
## Computation of distances for RFLP data
###############################################################################

## x: data.frame with RFLP data
## distfun: function to compute distance (cf. ?dist)
RFLPrefplot <- function(x, ref, distfun = dist, nrBands, mar.bottom = 5, 
                        cex.main = 1.2, cex.axis = 0.5, devNew = FALSE,
                        colBands, xlab = "", ylab = "molecular weight", 
                        ylim, ...){
    res <- RFLPdist2ref(x, ref, distfun, nrBands)
        
    x1 <- split(x, x$Sample)
    temp <- x1[names(x1) %in% rownames(res)]

    ref1 <- split(ref, ref$Sample)
    grp <- sapply(ref1, function(x) paste(x[1,"Taxonname"], " (", x[1,"Accession"], ")", sep = ""))
    ref.temp <- ref1[grp %in% colnames(res)]
    
    temp1 <- do.call("rbind", temp)
    ref.temp1 <- do.call("rbind", ref.temp)
    
    par(mar = c(mar.bottom, 4, 4, 2) + 0.1)
    if(missing(ylim)){
        lo <- 10*trunc(min(temp1$MW, ref.temp1$MW)/10)
        up <- 10*ceiling(max(temp1$MW, ref.temp1$MW)/10)
        if(up-lo < 1) up <- up + 1
        ylim <- c(lo, up)
    }
    At <- round(seq(ylim[1], ylim[2], length = 11), 0)
    
    grp1 <- paste(ref.temp1[, "Taxonname"], " (", ref.temp1[, "Accession"], ")", sep = "")
    if(devNew){
        par(mar = c(mar.bottom, 4, 4, 2) + 0.1)
    }else{
        par(mfrow = c(1,ncol(res)), mar = c(mar.bottom, 4, 4, 2) + 0.1)
    }
    for(i in 1:ncol(res)){
        if(devNew && i > 1){ 
            dev.new()
            par(mar = c(mar.bottom, 4, 4, 2) + 0.1)
        }
        temp.mw <- ref.temp1$MW[grp1 == colnames(res)[i]]
        reps <- length(temp.mw)
        samp.o <- order(res[,i])
        plot(NA, NA, xlim = c(0.25,length(temp)+2.75), ylim = c(lo, up),
            xlab = "", ylab = "molecular weight", axes = FALSE)
        axis(1, at = c(1,3:(length(temp)+2)), labels = c("Reference", rownames(res)[samp.o]), 
             las = 2, cex.axis = cex.axis)
        axis(2, at = At, labels = as.character(At), las = 2)
        abline(h = At, col = "grey", lty = 2)
        abline(h = temp.mw, col = "black", lty = 1)
        box()
        title(paste("Reference sample:", colnames(res)[i]), cex.main = cex.main)
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

        matlines(rbind(rep(0.75, reps), rep(1.25, reps)), 
                 rbind(temp.mw, temp.mw), 
                 lwd = 3, lty = 1, col = "black")
                    
        temp.o <- temp[samp.o]
        for(i in seq_along(temp.o)){
            temp.mw <- temp.o[[i]]$MW
            reps <- length(temp.mw)
            matlines(rbind(rep(i+2-0.25, reps), rep(i+2+0.25, reps)), 
                     rbind(temp.mw, temp.mw), 
                     lwd = 3, lty = 1, col = mycol[i])
        }
    }
    invisible()
}
