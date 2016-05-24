#' @export plot.pers
#' @title S3 Plotting Person - Item Map
#' @description S3 plotting method for object of class\code{"pers"}
#' @param x object of class\code{"pers"}
#' @param ra an integer, defining the (logit) range for y-axis
#' @param sortdif logical wether to order items by difficulty
#' @param main see \code{\link{plot}}
#' @param ylab see \code{\link{plot}}
#' @param itemNames logical wether to use itemnames in the resulting plot
#' @param fillCol color for bar filling of the ability histogram
#' @param lineCol color for bar lines of the ability histogram
#' @param cex see \code{\link{text}}
#' @param pos see \code{\link{text}}
#' @param breaks see \code{\link{hist}}
#' @param pch see \code{\link{points}}
#' @param ... other parameters passed to \code{\link{hist}} and \code{\link{text}}.

plot.pers<-function(x, ra=NULL, sortdif=FALSE, main=NULL, ylab="Logits", itemNames=TRUE,fillCol="grey60", lineCol="grey40", cex=.7, pos=4, breaks = "Sturges", pch = 1, ...){ 
  if(length(main)==0){
    main<-paste("Person - Item Map ","\n ", deparse(substitute(x)) , sep="")
    #main<-deparse(substitute(x))
  }
  pers_obj <- x
if(sortdif==TRUE){
  threshold <- pers_obj$pair$threshold
  sigma <- pers_obj$pair$sigma
  spaltennamen <- colnames(threshold)
  #####
  threshold <- as.matrix(threshold[order(sigma), ])
  colnames(threshold) <- spaltennamen
  sigma <- sort(sigma)
  cat("(ordered by location) \n")
}else{
  threshold <- pers_obj$pair$threshold
  sigma <- pers_obj$pair$sigma}
  

#---------------------------------
hist_obj <- hist(pers_obj$pers$WLE, plot = FALSE, breaks=breaks)#  , ...

binWidth <- hist_obj$breaks[2] - hist_obj$breaks[1]
xscale <- 1 * 0.90 / max(hist_obj$density)
xpos <- max(hist_obj$density)*1.2

def_par <- par(no.readonly = TRUE) # save default, for resetting...

pw <- layout(matrix(c(1,2),1,2), widths = c(2,3), TRUE) # layout.show(pw)

par(oma=c(0,0,4,0))
par(mar=c(3,4,0,0))
#-------------
n <- length(hist_obj$density)
x.l <- rep(xpos, n)*xscale
x.r <- x.l-(hist_obj$density * xscale)
y.b <- hist_obj$breaks[1:n]
y.t <- hist_obj$breaks[2:(n + 1)]
#-------------
if(length(ra)==0){ra <- max(abs(hist_obj$breaks))}

plot(x=0:(xpos*xscale), y=c(-ra,ra), type="n", bty="n", xaxt="n", ylab=ylab)
# plot(x=0:(xpos*xscale), y=range(hist_obj$breaks), type="n", bty="n", xaxt="n", ylab=ylab)
rect(xleft = x.l, ybottom = y.b, xright = x.r, ytop = y.t, col = fillCol, border = lineCol)
#---------------------------'
par(mar=c(3,1,0,3))
#-------------
# plot(x=1:(length(sigma)+1),y=seq(min(hist_obj$breaks),max(hist_obj$breaks),length.out=(length(sigma)+1)),type="n",bty="n",ylab="" ,yaxt="n",xaxt="n" ,xlab="items")
plot(x=1:(length(sigma)+1),y=seq(-ra,ra,length.out=(length(sigma)+1)),type="n",bty="n",ylab="" ,yaxt="n",xaxt="n" ,xlab="items")

if (itemNames==TRUE){text(x=1:length(sigma),y=sigma, labels=names(sigma),cex=cex,pos=pos, ...)}

#-------------
for ( i in 1:dim(threshold)[1] ){
  lines(x=rep(i,length(na.omit(threshold[i,])) ), y= na.omit(threshold[i,]), type = "l")
  points(x=rep(i,length(na.omit(threshold[i,])) ), y= na.omit(threshold[i,]), pch=pch)
}
#-------------
title(main = main, sub = NULL, xlab = NULL, ylab = NULL,line = 1, outer = TRUE)

par(def_par) # resetting parameter

}
