plot.steeptest <- function (x,...) {
  
  if ((is.null(x$names)) && (length(x$names)<=15)) names <- paste('Ind.',1:nrow(x$matdom))
  if ((is.null(x$names)) && (length(x$names)>15)) names <- 1:nrow(x$matdom)
  if ((nchar(x$names)<10) && (length(x$names)<=15)) names <- x$names
  if ((nchar(x$names)>10) || (length(x$names)>15)) names <- 1:nrow(x$matdom)
  
  SortNormDS <- sort(x$NormDS,decreasing=TRUE,index.return=TRUE)
  names <- names[SortNormDS$ix]
  NormDS <- array(SortNormDS$x,dim=nrow(x$matdom),dimnames=list(names))
  rnk <- array(1:nrow(x$matdom))
 
  if (exists('xlab')==FALSE) xlab <- c("Individuals in Rank Order")
  if (exists('ylab')==FALSE) ylab <- c("Normalized David's Scores")
  if (exists('main')==FALSE) main <- paste("NormDS (based on",  if (x$method == 'Dij')
    "Dij)" else "Pij)", "plotted against rank order") 
  
  plot(rnk,NormDS, type = "b", lty = 2, ylim = c(0,(nrow(x$matdom)-1)), xlab=xlab, ylab=ylab,
       main=main, pch=15, col="blue",axes=FALSE,xaxs = "r",yaxs = "i")

  axis(1, at=1:length(names),names,tcl = 0.3, adj=c(0.5))
  axis(2, at=0:(nrow(x$matdom)-1),tcl = 0.3, srt = 1)
  box()
  myline.fit <- lm(NormDS ~ rnk)

  abline(myline.fit,col="green",lwd=2)
  lgnd <- paste("Fitted line:  Y = ",format(myline.fit$coefficients[2],digits=4),
  "X", "+" ,format(myline.fit$coefficients[1],digits=4))
  legend("topright",as.graphicsAnnot(lgnd),lty=1,col="green",lwd=2,box.lty=0)
  invisible()  
}
