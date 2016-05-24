plots.confints.bootpls = function (ic_bootobject, indices = NULL, legendpos = "topleft", 
    prednames = TRUE, articlestyle = TRUE, xaxisticks=TRUE, ltyIC=c(2, 4, 5, 1), colIC=c("darkgreen", "blue", "red", "black"), typeIC, las=par("las"), mar, mgp, ...) 
{  
    if(missing(typeIC)){
      if(attr(ic_bootobject, "typeBCa")){
        typeIC <- c("Normal", "Basic", "Percentile", "BCa")
      }
      else {
        typeIC <- c("Normal", "Basic", "Percentile")
      }
    }
    if((!attr(ic_bootobject, "typeBCa"))&("BCa" %in% typeIC)){stop("BCa intervals were not computed, hence cannot be plotted.")}
    if(length(ltyIC)<length(typeIC)){ltyIC <- rep_len(ltyIC,length(typeIC))}
    if(length(colIC)<length(typeIC)){colIC <- rep_len(colIC,length(typeIC))}
    nr <- nrow(ic_bootobject)
    if (is.null(indices)) {
        indices <- 1:nr
    }
    plotpos <- (1:nr)[1:length(indices)]
    if (articlestyle) {
      oldparmar <- par("mar")
      oldparmgp <- par("mgp")
      if(missing(mar)){mar=c(2, 2, 1, 1) + 0.1}
        if(missing(mgp)){mgp=c(2, 1, 0)}
        par(mar = mar); par(mgp = mgp)
    }
    plot(c(1, 1), xlab = "", ylab = "", type = "n", xlim = c(1, 
        length(indices) + 0.5), ylim = c(min(ic_bootobject[indices,]), 
        max(ic_bootobject[indices,])), xaxt = "n", ...)
        legendtxt <- NULL
        indictypeIC <- rep(FALSE,4)
        nbIC <- 0
        if ("Normal" %in% typeIC){
        indictypeIC[1] <- TRUE
        arrows(plotpos + nbIC*0.15, ic_bootobject[indices, 1], plotpos, ic_bootobject[indices, 
            2], lend = "butt", lwd = 2, lty = ltyIC[1], col = colIC[1], 
            code = 3, angle = 90, length = 0.1)
            legendtxt <- c(legendtxt,"Normal")
        nbIC <- nbIC+1
            }
        if ("Basic" %in% typeIC){        
        indictypeIC[2] <- TRUE
        arrows(plotpos + nbIC*0.15, ic_bootobject[indices, 3], plotpos + 
            nbIC*0.15, ic_bootobject[indices, 4], lend = "butt", lwd = 2, 
            lty = ltyIC[2], col = colIC[2], code = 3, angle = 90, length = 0.1)
            legendtxt <- c(legendtxt,"Basic")
        nbIC <- nbIC+1
            }
        if ("Percentile" %in% typeIC){        
        indictypeIC[3] <- TRUE
        arrows(plotpos + nbIC*0.15, ic_bootobject[indices, 5], plotpos + 
            nbIC*0.15, ic_bootobject[indices, 6], lend = "butt", lwd = 2, 
            lty = ltyIC[3], col = colIC[3], code = 3, angle = 90, length = 0.1)
            legendtxt <- c(legendtxt,"Percentile")
        nbIC <- nbIC+1
        }
        if (("BCa" %in% typeIC)&(attr(ic_bootobject, "typeBCa"))){
        indictypeIC[4] <- TRUE
        arrows(plotpos + nbIC*0.15, ic_bootobject[indices, 7], plotpos + 
            nbIC*0.15, ic_bootobject[indices, 8], lend = "butt", lwd = 2, 
            lty = ltyIC[4], col = colIC[4], code = 3, angle = 90, length = 0.1)
            legendtxt <- c(legendtxt,"BCa")     
        nbIC <- nbIC+1
        }
        if (prednames) {
            if(xaxisticks){
              axis(1, at = plotpos + (nbIC-1)*0.15/2, labels = rownames(ic_bootobject)[indices], las=las)
            }
            else
            {
              axis(1, at = plotpos + (nbIC-1)*0.15/2, labels = rownames(ic_bootobject)[indices],lwd.ticks=0, las=las)
            }
        }
        else {
            if(xaxisticks){
              axis(1, at = plotpos + (nbIC-1)*0.15/2, labels = paste("x", (1:nr)[indices], sep = ""), las=las)
            }
            else
            {
              axis(1, at = plotpos + (nbIC-1)*0.15/2, labels = paste("x", (1:nr)[indices], sep = ""),lwd.ticks=0, las=las)
            }
        }
        abline(h = 0, lty = 3, lwd = 2)
        legend(legendpos, legend = legendtxt, lty = ltyIC[indictypeIC], col = colIC[indictypeIC], lwd = 2)
    if (articlestyle) {
      par(mar=oldparmar)
      par(mgp=oldparmgp)
    }
}
