##' Plot ellipses for two-dimensional data
##' 
##' The function plots ellipses for different categories from two-dimensional
##' data.
##' 
##' 
##' @param x A two-columned matrix of data
##' @param labs An optional vector of labels, parallel to 'data'
##' @param chars An optional vector of labels, parallel to 'data'. If this
##' argument is specified these labels will be plotted rather than the labels
##' in 'labs'.
##' @param formant If TRUE) then the data is negated and the axes are switched
##' so that, for formant data, the plot is made with decreasing F2 on the
##' x-axis and decreasing F1 on the y-axis.
##' @param scaling Either "mel" or "bark" for mel or bark scaling of the data
##' @param prob A single numeric vector greater than zero and less than 1
##' representing the confidence interval of the ellipse contours. Defaults to
##' 0.95
##' @param nsdev Defines the length of the major and minor axes of the ellipses
##' in terms of the standard deviation of the data and overrides the prob
##' argument.
##' @param dopoints If TRUE) character labels (from 'labs' or 'chars') are
##' plotted for each data point
##' @param doellipse If TRUE, ellipses are drawn on the plot. If FALSE, no
##' ellipses are drawn and, if 'dopoints' is also FALSE, 'centroids' is set to
##' T
##' @param centroid One label for each ellipse is drawn
##' @param axes If TRUE axes are drawn on the plot.
##' @param xlim A vector of two numeric values giving the range of the x-axis.
##' @param ylim A vector of two numeric values giving the range of the y-axis.
##' @param col If colour is TRUE) the ellipses and labels will be plotted in
##' different colours
##' @param lty If linetype is TRUE) the ellipses will be plotted with different
##' linetypes.  This is useful for plots that will be printed.
##' @param lwd A code passed to the lwd argument in plotting functions. 'lwd'
##' can be either a single element numeric vector, or its length must be equal
##' to the number of unique types in labs. For example, if lwd=3 and if labs =
##' c("a", "b", "a", "c"), then the output is c(3, 3, 3, 3). Alternatively, if
##' lwd = c(2,3,1), then the output is c(2, 3, 2, 1) for the same example. The
##' default is NULL in which case all lines are drawn with lwd=1
##' @param ... graphical options \link{par}
##' @return NULL
##' @author Jonathan Harrington, Steve Cassidy
##' @seealso \code{\link{dcut}}
##' @keywords dplot
##' @examples
##' 
##' 
##'    data(vowlax)
##'    data <- cbind(vowlax.df$F1,vowlax.df$F2)
##'    phonetic = vowlax.l
##'    word = vowlax.word
##' 
##'    eplot(data, phonetic)
##'     
##'   
##'    eplot(data, phonetic, form=TRUE, main="F1 x F2 plane", centroid=TRUE)
##'    eplot(data, phonetic, form=TRUE, main="F1 x F2 plane", dopoints=TRUE)
##'    eplot(data, phonetic, form=TRUE, main="F1 x F2 plane in Bark", 
##'          dopoints=TRUE, scaling="bark")
##'    eplot(data, phonetic, form=TRUE, main="F1 x F2 plane in Bark b/w with linetype", 
##'          col=FALSE, lty=TRUE, dopoints=TRUE, scaling="bark") 
##'    eplot(data, phonetic, form=TRUE, main="F1 x F2 plane", 
##'          doellipse=FALSE, dopoints=TRUE)
##'    eplot(data, phonetic, form=TRUE, dopoints=TRUE, 
##'          prob=0.5, main="F1 x F2 plane, 50% confidence intervals")
##'    eplot(data, phonetic, form=TRUE, dopoints=TRUE, 
##'          nsdev=2, main="F1 x F2 plane, 2 standard deviations")
##'    
##' 
##'    temp <- phonetic %in% c("a", "O")
##'    eplot(data[temp,], phonetic[temp], form=TRUE,  main="F1 x F2 [A] and [O] only", centroid=TRUE)
##'    
##' 
##'    temp <- phonetic=="O"
##'    eplot(data[temp,], phonetic[temp], word[temp], form=TRUE, 
##'          dopoints=TRUE, main="[O] only showing word labels")  
##'       
##'    
##' 
##' 
##' 
##' @export eplot
`eplot` <- function (x, labs, chars, formant = FALSE, scaling = "linear", 
                     prob = 0.95, nsdev = NULL, dopoints = FALSE, 
                     doellipse = TRUE, centroid = FALSE,  axes = TRUE, 
                     xlim, ylim, col = TRUE, lty = FALSE,  lwd = NULL, ...) 
{
  ocall <- match.call()
  if (is.null(nsdev)) 
    nsdev <- sqrt(qchisq(prob, 2))
  
  if (missing(labs)) 
    labs <- rep(".", nrow(x))
  if (!doellipse & !dopoints) 
    centroid <- TRUE
  if (nrow(x) != length(labs)) 
    stop("x and labels don't match")
  if (ncol(x) != 2) 
    stop("Eplot needs 2 dimensional x")
  if (!missing(chars)) 
    if (length(labs) != length(chars)) 
      stop("Length of chars must match that of labs")
  if (scaling == "mel") 
    x <- mel(x)
  if (scaling == "bark") 
    x <- bark(x)
  if (formant) {
    x <- cbind(-x[, 2], -x[, 1])
    if (!missing(xlim)) 
      xlim <- -rev(xlim)
    if (!missing(ylim)) 
      ylim <- -rev(ylim)
  }
  col.lty <- mu.colour(labs, col, lty, lwd)
  lty <- col.lty$linetype
  linewidth <- col.lty$lwd
  uniqlabels <- unique(labs)
  emat <- nums <- cen <- k <- l <- NULL
  for (j in uniqlabels) {
    temp <- labs == j
    mat <- x[temp, , drop = FALSE]
    if (nrow(mat) > 2) {
      evals <- eigen(var(mat))
      m1 <- mean(mat[, 1])
      m2 <- mean(mat[, 2])
      e <- ellipse(m1, m2, sqrt(evals$values[1]) * nsdev, 
                   sqrt(evals$values[2]) * nsdev, aperm(evals$vectors, 
                                                        c(2, 1)))
    }
    else {
      cat("Too few x points for label ", j, " will plot a point or a line\n")
      m1 <- mean(mat[, 1])
      m2 <- mean(mat[, 2])
      e <- mat
    }
    nums <- c(nums, nrow(e))
    emat <- rbind(emat, e)
    k <- c(k, col.lty$legend$col[match(j, col.lty$legend$lab)])
    l <- c(l, col.lty$legend$lty[match(j, col.lty$legend$lab)])
    linewidth <- c(linewidth, col.lty$legend$lwd[match(j, 
                                                       col.lty$legend$lab)])
    if (centroid) 
      cen <- rbind(cen, cbind(m1, m2))
  }
  if (doellipse) {
    if (missing(xlim)) 
      xlim <- range(c(emat[, 1], x[, 1]))
    if (missing(ylim)) 
      ylim <- range(c(emat[, 2], x[, 2]))
  }
  else {
    if (missing(xlim)) 
      xlim <- range(x[, 1])
    if (missing(ylim)) 
      ylim <- range(x[, 2])
  }
  rightlim <- cumsum(nums)
  leftlim <- cumsum(nums) - (nums - 1)
  rowmarker <- cbind(leftlim, rightlim)
  for (j in 1:nrow(rowmarker)) {
    lowerlim <- rowmarker[j, 1]
    upperlim <- rowmarker[j, 2]
    if (doellipse) {
      graphics::plot(emat[lowerlim:upperlim, ], type = "l", axes = FALSE, 
           xlim = xlim, ylim = ylim, col = k[j], 
           lty = as.numeric(l[j]), lwd = as.numeric(linewidth[j]), xlab="", ylab="", main="")
    }
    else {
      graphics::plot(emat[lowerlim:upperlim, ], type = "n", axes = FALSE, 
           xlim = xlim, ylim = ylim, col = k[j], 
           lty = as.numeric(l[j]), lwd = as.numeric(linewidth[j]), xlab="", ylab="", main="")
    }
    if (dopoints) {
      centroid <- FALSE
      singlelab <- uniqlabels[j]
      temp <- labs == singlelab
      if (!missing(chars)) 
      {
        if(is.numeric(chars))
          graphics::points(x[temp, 1], x[temp, 2], pch=chars[temp], 
                 col = k[j])
        else
          graphics::text(x[temp, 1], x[temp, 2], chars[temp], 
               col = k[j])
      }
      else graphics::text(x[temp, 1], x[temp, 2], labs[temp], 
                col = k[j])
    }
    if (centroid) {
      singlelab <- uniqlabels[j]
      graphics::text(cen[j, 1], cen[j, 2], singlelab, col = k[j])
    }
    if (j < nrow(rowmarker)) 
      graphics::par(new = TRUE)
  }
  graphics::par(col = 1)
  
  if (axes) {
    if (formant) {
      xaxp <- graphics::par("xaxp")
      yaxp <- graphics::par("yaxp")
      xat <- seq(xaxp[1], xaxp[2], length.out = xaxp[3] + 
                   1)
      yat <- seq(yaxp[1], yaxp[2], length.out = yaxp[3] + 
                   1)
      graphics::axis(1, at = xat, labels = -xat)
      graphics::axis(2, at = yat, labels = -yat, srt = 90)
    }
    else {
      graphics::axis(1)
      graphics::axis(2)
    }
  }
  graphics::title(...)
  graphics::box(...)
}











##' Calculate ellipse coordinates
##' 
##' Calculates ellipse coordinates for eplot
##' 
##' 
##' @param x X coordinate of center
##' @param y y coordinate of center
##' @param rx Radius in the x direction
##' @param ry Radius in the y direction
##' @param orient Orientation, in radians. The angle of the major axis to the x
##' axis.
##' @param incr The increment between points, in degrees.
##' @return A matrix of x and y coordinates for the ellipse.
##' @seealso \code{\link{eplot}}
##' @keywords misc
##' @export ellipse
"ellipse"<- function(x, y, rx, ry, orient, incr = 360/100)
{
  rincr <- radians(incr)
  theta <- seq(0, 2 * pi, rincr)
  xcoord <- rx * cos(theta)
  ycoord <- ry * sin(theta)
  mat <- cbind(xcoord, ycoord)
  mat <- mat %*% orient
  mat[, 1] <- mat[, 1] + x
  mat[, 2] <- mat[, 2] + y
  mat
}









##' polygonplot
##' 
##' plots a polygon
##' 
##' 
##' @param data data matrix
##' @param labels labels
##' @param order order
##' @param formant formant TRUE or FALSE transposes the axes
##' @param axes axes
##' @param xlab xlab
##' @param ylab ylab
##' @param main main
##' @param xlim xlim
##' @param ylim ylim
##' @keywords internal
##' @export polygonplot
"polygonplot" <- function(data, labels, order,
                          formant=TRUE, axes=TRUE,
                          xlab="", ylab="",
                          main = "", xlim, ylim)
{
  
  if( ncol(data) > 2 ) {
    data <- data[,1:2]
  }
  if( ncol(data) != 2 ) {
    stop( "polygonplot() requires two columns of data" )
  }
  
  if(formant)
    data <- cbind(-data[, 2], -data[, 1])
  
  
  points <- NULL
  for( l in order ) {
    tmp <- matrix(data[labels==l],ncol=2)
    points <- rbind( points, apply(tmp, 2, mean) )
  }
  
  graphics::plot( points, type="b", pch=" ", axes=FALSE, xlab="", ylab="" )
  graphics::text( points, order, axes=FALSE, , xlab="", ylab="" )
  
  graphics::par(col = 1)
  graphics::box()
  if(axes) {
    if(formant) {
      if(missing(xlab))
        xlab <- "F2"
      if(missing(ylab))
        ylab <- "F1"
      xaxp <- graphics::par("xaxp")
      yaxp <- graphics::par("yaxp")
      xat <- seq(xaxp[1], xaxp[2], length.out = xaxp[3] + 1)
      yat <- seq(yaxp[1], yaxp[2], length.out = yaxp[3] + 1)
      graphics::axis(1, at = xat, labels =  - xat)
      graphics::axis(2, at = yat, labels =  - yat, srt = 90)
    }
    else {
      graphics::axis(1)
      graphics::axis(2)
    }
  }
  graphics::title(main = main, xlab = xlab, ylab = ylab)
}
