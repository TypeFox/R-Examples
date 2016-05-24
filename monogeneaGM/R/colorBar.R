#' Add a color bar
#'
#' This function produces a color bar. 
#' @param colpalette a color palette supplied by the user; defaults to the red-white-blue palette if left undefined
#' @param min the smallest numeric value corresponding to the color on the extreme left of the color palette
#' @param max the largest numeric value corresponding to the color on the extreme right of the color palette
#' @param nticks number of tick lines on the color bar
#' @param ticks spacing for the tick marks on the color bar
#' @param tit title for the color bar 
#' @details This R code is a based on John Colby's (2011) \code{color.bar} code. A user-defined color palette can be generated using the 
#' \code{colorRampPalette} function. See www.colorbrewer2.org for interesting color palette options.
#' @seealso \code{\link{pcloadhm}}, \code{\link{colorRampPalette}}
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Colby J. (2011). Color bar legends for neuroimaging in R. Available at http://www.colbyimaging.com/wiki/statistics/color-bars.
#'
#' Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#' @examples
#' #define a color scale
#' tones <- c("#99000D", "#FB6A4A", "white", "#6BAED6", "#084594")
#' multihue <- colorRampPalette(tones, space = "Lab")(101)
#'
#' #create a two-panel figure with width ratio of 4:1
#' nf <- layout(matrix(c(1,1,1,1,2),2,5,byrow=TRUE))
#' layout.show(nf)
#'
#' #mapping a matrix of randomly chosen numbers between -1 and 1 to
#' #colors in the color scale 
#' h <- matrix(runif(100, -1, 1),10,10)
#' image(h,col=multihue,xaxt="n", yaxt="n")
#'
#' #add color bar
#' par(mar=c(5,5,5,2))
#' colorBar(colpalette=multihue, min=-1,max=1)
#'

 colorBar <-function(colpalette=NULL, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), 
 tit="") {
    
	 if(is.null(colpalette)) {
      scale <- (length(colpalette)-1)/(max-min)
	 rwb <- c("#99000D","#FB6A4A","white","#6BAED6","#084594")
	 colpalette<- colorRampPalette(rwb, space="Lab")(101)
	 }

    scale <- (length(colpalette)-1)/(max-min)
    plot(c(0,10), c(min,max), type="n", bty="n", xaxt="n", xlab="", yaxt="n", ylab="", main=tit)
    axis(2, ticks, las=1)
    for (i in 1:(length(colpalette)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=colpalette[i], border=NA)
    }
 }
