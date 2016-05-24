plotdistmat <-
function(distmat, colvec, coltext, pos="topleft", labels=NULL, numbers=TRUE, ...) {
  tot.inf <- nrow(distmat)
  maxD <- max(distmat)
  col.length <- length(colvec)
  if (length(labels)==0) {
    labels <- 1:tot.inf
  }
  plot(NULL, xlim=c(0,tot.inf), ylim=c(0,tot.inf), main="", xlab="", 
       ylab="", bty="n", xaxs="i", yaxs="i", xaxt="n", yaxt="n")
  if (pos=="topleft") {
    mtext("Genetic distance", line=2)
    axis(2, at=(1:tot.inf)-0.5, labels=rev(labels), las=1, cex.axis=0.5)
    axis(3, at=(1:tot.inf)-0.5, labels=rev(labels), las=2, cex.axis=0.5)
    for (i in 1:tot.inf) {
      for (j in 1:(tot.inf-i+1)) {
        colz <- colvec[floor((col.length-1)*(distmat[i,tot.inf-j+1])/maxD)+1]
        textcol <- coltext[floor((col.length-1)*(distmat[i,tot.inf-j+1])/maxD)+1]
        rect(j-1, tot.inf-i, j, tot.inf-i+1, border=NA, col=colz)
        if (numbers) {
          text(j-0.5, tot.inf-i+0.5, distmat[i,tot.inf-j+1], adj=c(0.5,NA), col=textcol, cex=0.5)
        }
      }
    }
  } else if (pos=="topright") {
    mtext("Genetic distance", line=2)
    axis(3, at=(1:tot.inf)-0.5, labels=labels, las=1, cex.axis=0.5)
    axis(4, at=(1:tot.inf)-0.5, labels=rev(labels), las=2, cex.axis=0.5)
    for (i in 1:tot.inf) {
      for (j in i:tot.inf) {
        colz <- colvec[floor((col.length-1)*(distmat[i,j])/maxD)+1]
        textcol <- coltext[floor((col.length-1)*(distmat[i,j])/maxD)+1]
        rect(j-1, tot.inf-i, j, tot.inf-i+1, border=NA, col=colz)
        if (numbers) {
          text(j-0.5, tot.inf-i+0.5, distmat[i,j], adj=c(0.5,NA), col=textcol, cex=0.5)
        }
      }
    }
  } else if (pos=="bottomleft") {
    axis(2, at=(1:tot.inf)-0.5, labels=rev(labels), las=1, cex.axis=0.5)
    axis(1, at=(1:tot.inf)-0.5, labels=labels, las=1, cex.axis=0.5)        
    for (j in 1:tot.inf) {
      for (i in j:tot.inf) {
        colz <- colvec[floor((col.length-1)*(distmat[i,j])/maxD)+1]
        textcol <- coltext[floor((col.length-1)*(distmat[i,j])/maxD)+1]
        rect(j-1, tot.inf-i, j, tot.inf-i+1, border=NA, col=colz)
        if (numbers) {
          text(j-0.5, tot.inf-i+0.5, distmat[i,j], adj=c(0.5,NA), col=textcol, cex=0.5)
        }
      }
    }
  } else if (pos=="bottomright") {
    axis(4, at=(1:tot.inf)-0.5, labels=rev(labels), las=1, cex.axis=0.5)
    axis(1, at=(1:tot.inf)-0.5, labels=rev(labels), las=1, cex.axis=0.5)  
    for (i in 1:tot.inf) {
      for (j in (tot.inf-i+1):tot.inf) {
        colz <- colvec[floor((col.length-1)*(distmat[i,tot.inf-j+1])/maxD)+1]
        textcol <- coltext[floor((col.length-1)*(distmat[i,tot.inf-j+1])/maxD)+1]
        rect(j-1, tot.inf-i, j, tot.inf-i+1, border=NA, col=colz)
        if (numbers) {
          text(j-0.5, tot.inf-i+0.5, distmat[i,tot.inf-j+1], adj=c(0.5,NA), col=textcol, cex=0.5)
        }
      }
    }
  }
}
