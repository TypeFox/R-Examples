#----- Function for plotting matrices ----- #
# Curtesy of Chris Seidel
# Available at http://www.phaget4.org/R/image_matrix.html
# optional arguments: myImagePlot(m, xLabels, yLabels, zlim, title=c("my title"))

.plotMatrix <- function(x, correl = FALSE, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  if( is.null(Lst$ylab) ){
    ylab <- ""
  } else {ylab <- c(Lst$ylab)}
  if( is.null(Lst$xlab) ){
    xlab <- ""
  } else {xlab <- c(Lst$xlab)}
  if( is.null(Lst$scaleLab) ){
    scaleLab <- ""
  } else {scaleLab <- c(Lst$scaleLab)}
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  if(correl){
    ColorRamp <- colorRampPalette(c("red", "salmon", "white", "royalblue", "navy"))(256)
    min <- -1
    max <- 1
  } else{
    ColorRamp <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(256)
  }
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  if(correl) {
    diag(x) <- 0
    x[upper.tri(x)] <- 0
  }
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(4.5, 4.5, 2.5, 2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab = xlab, ylab = ylab, axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="", main = scaleLab,
        xaxt="n")
  
  layout(1)
}