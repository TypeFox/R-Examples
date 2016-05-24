### plot.riskReclassification.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Sep 24 2015 (19:26) 
## Version: 
## last-updated: Oct  3 2015 (16:26) 
##           By: Thomas Alexander Gerds
##     Update #: 7
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' @export
plot.riskReclassification <- function(x,xlim=c(0,100),ylim=c(0,100),xlab,ylab,grid=TRUE,grid.col=gray(0.9),...){
    if (missing(xlab)) xlab <- paste("Risk (%):",names(dimnames(x$reclassification))[[1]])
    if (missing(ylab)) ylab <- paste("Risk (%):",names(dimnames(x$reclassification))[[2]])
    plot(x$predictedRisk[[1]],
         x$predictedRisk[[2]],
         axes=FALSE,
         xlim=xlim,
         ylim=ylim,
         xlab=xlab,
         ylab=ylab,
         ...)
    axis(1,at=x$cuts)
    axis(2,at=x$cuts)
    if (grid==TRUE)
        abline(h = x$cuts, v = x$cuts, col = gray(0.9))
}


#----------------------------------------------------------------------
### plot.riskReclassification.R ends here
