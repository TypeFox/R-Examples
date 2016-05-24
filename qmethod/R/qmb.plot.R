qmb.plot <- function(qmbsum, type=c("zsc", "loa"), nfactors, cex = 0.7, cex.leg=0.8, errbar.col= "black", lwd=1, lty=1, vertdist = 0.2, limits=NULL, r.names=NA, sort=c("none", "difference", "sd"), sbset=NULL, leg.pos="topleft", bty = "n", plot.std = TRUE, pch= NULL, col=NULL, grid.col="gray", ...) {
  if(type == "loa") {
    boloa <- qmbsum[[1]]
    db <- boloa[ ,c(grep("loa", names(boloa)), grep("SE", names(boloa)), grep("std", names(boloa)))]
    item <- "Q-sort"
    values <- "Factor loading"
    if(is.null(limits)) limits <- c(-1.0,1.0)
  }
  
  if(type == "zsc") {
    boloa <- qmbsum[[2]]
    db <- boloa[ ,c(grep("zsc.bts", names(boloa)), grep("SE", names(boloa)),   grep("std", names(boloa)))]
    item <- "Statement"
    values <- "z-score"
    if(is.null(limits)) {
      zscs <- grep("zsc.bts", names(db))
      SEs <- grep("SE", names(db))
      lms.down <- db[,zscs] - db[,SEs]
      lms.up <- db[,zscs] + db[,SEs]
      limits <- c(floor(min(lms.down)), ceiling(max(lms.up)))
    }
  }
  if(is.numeric(sbset)) db <- db[c(1:min(nrow(db), sbset)), ]
  nitems <- nrow(db)
  if(length(r.names) == nrow(db)) rownames(db) <- r.names
  
  if(sort == "sd") {
    sds <- apply(db[,(1+nfactors):(2*nfactors)], 1, sum)
    db <- db[order(sds), ]
  }
  if(sort == "difference") {
    sds <- abs(apply(db[,(1:nfactors)], 1, sd))
    db <- db[order(sds), ]
  }
  
  #Plotting parameters
  db$position <- c(1:nitems)
  if(is.null(col)) {
    colegend=c(rep("black", nfactors), rep("white", 3))
    dot.col="black"
    } else {
      colegend=c(col[1:nfactors], rep("white", 3))
      dot.col=col[1:nfactors]
    }
  if(is.null(pch)) pich=array(c(20, 15, 17, 18)) else pich=pch[1:nfactors]
  if(is.null(pch)) {
    pitx=array(c(21, 22, 24, 23)) 
    } else if(plot.std) {
      if (length(pch) >= 2*nfactors) {
        pitx=pch[(nfactors+1):(2*nfactors+1)]
        } else stop("The vector of symbols provided in 'pch' needs to be at least twice the length of the number of factors, in order to contain (a) a set of symbols for the bootstrap values and (b) a different set of symbols for the standard values.")
    }
  i=1
  # Plot:
  dotchart(db[,i], labels = rownames(db), pch=pich[i], 
           xlim=limits,  
           xlab=values, lcolor="white",
           lwd = lwd, cex=cex, color=dot.col[i], ...)
  mtext(item, side=2, line=1.5, cex=cex, ...)
  # Error bars:
  segments(x0=db[,i], y0=db[,"position"], 
           x1=db[,i]+db[,nfactors+i], 
           y1=db[,"position"], lwd = lwd, lty = lty, col = errbar.col, cex=cex, ...)
  segments(x0=db[,i], y0=db[,"position"], 
           x1=db[,i]-db[,nfactors+i], 
           y1=db[,"position"], lwd = lwd, lty = lty, col = errbar.col, cex=cex, ...)
  # Replot points, for them to be on top of error bars:
  points(x=db[,i], db[,"position"]+(vertdist*(i-1)), 
         pch = pich[i], type = "p", lwd = lwd, 
         cex=cex, col=dot.col[i], ...)
  if(plot.std) {
    points(x=db[,(2*nfactors)+i], db[,"position"], pch = pitx[i], 
           type = "p", lwd = lwd, cex=cex, col=dot.col[i], ...)
  }
  # Plot 2nd and subsequent factors
  for (i in 2:nfactors) {
    # Error bars:
    segments(x0=db[,i], y0=db[,"position"]+(vertdist*(i-1)), 
             x1=db[,i]+db[,nfactors+i], 
             y1=db[,"position"]+(vertdist*(i-1)), lwd = lwd, 
             lty = lty, col = errbar.col, cex=cex, ...)
    segments(x0=db[,i], y0=db[,"position"]+(vertdist*(i-1)), 
             x1=db[,i]-db[,nfactors+i], 
             y1=db[,"position"]+(vertdist*(i-1)), lwd = lwd, 
             lty = lty, col = errbar.col, cex=cex, ...)
    # Points:
    points(x=db[,i], db[,"position"]+(vertdist*(i-1)), 
           pch = pich[i], type = "p", lwd = lwd, 
           cex=cex, col=dot.col[i], ...)
    if(plot.std) {
      points(x=db[,(2*nfactors)+i], 
             db[,"position"]+(vertdist*(i-1)), 
             pch = pitx[i], type = "p", lwd = lwd, 
             cex=cex, col=dot.col[i], ...)
    }
  }
  abline(v=seq(floor(limits[1]), ceiling(limits[2]), 0.5), col=grid.col, lty="dotted", lwd = lwd, ...)
  abline(h=c(0.7:(nitems+0.7)), col=grid.col, lty="dotted", lwd = lwd, ...)
  if(plot.std) leg.length <- 1:(nfactors+2) else leg.length <- 1:nfactors
  legend(leg.pos, legend=c(paste0("Factor ", 1:nfactors), "Empty symbol: standard", "Filled symbol: bootstrap")[leg.length], pch=c(pich[1:nfactors], 0, 0)[leg.length], cex=cex.leg*cex, pt.cex=cex, col=colegend, bg=NA, bty=bty, ...)
}