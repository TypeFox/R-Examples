##' @export 
plot.Cindex <- function(x,ylim=c(.4,1),xlim=c(0,x$maxtime),abline=TRUE,xlab="Time",ylab="Concordance index",...){
  argList <- match.call(expand.dots=TRUE)
  argList[[1]] <- as.name("list")
  argList <- eval(argList,parent.frame())
  argList <- c(list("what"=switch(x$splitMethod$internal.name,
                        "noPlan"={"AppCindex"},
                        paste(x$splitMethod$internal.name,"Cindex",sep=""),
                        xlab=xlab,ylab=ylab)),argList)
  argList$ylim <- ylim
  argList$xlim <- xlim
  argList$ylab <- ylab
  argList$xlab <- xlab
  argList$x$exact <- FALSE
  do.call("plot.pec", argList)
  if (abline==TRUE)
    abline(h=.5,col="gray",lty=3,lwd=3,xpd=FALSE)
}
