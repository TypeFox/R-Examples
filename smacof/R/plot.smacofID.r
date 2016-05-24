# plot method for all smacof objects

plot.smacofID <- function(x, plot.type = "confplot", plot.dim = c(1,2), bubscale = 1, col = 1, 
                          label.conf = list(label = TRUE, pos = 3, col = 1), identify = FALSE, 
                          type = "p", pch = 20,  asp = 1, main, xlab, ylab, xlim, ylim, ...)

# x ... object of class smacofID
# plot.type ... types available: "confplot", "bubbleplot", "stressplot"
# Shepard plot and resplot are performed over sum of distances
  
{
  ## --- check label lists
  if (is.null(label.conf$label)) label.conf$label <- TRUE
  if (is.null(label.conf$pos)) label.conf$pos <- 3
  if (is.null(label.conf$col)) label.conf$col <- 1
  if (is.null(label.conf$cex)) label.conf$cex <- 0.8
  if (identify) label.conf$label <- FALSE
  
  x1 <- plot.dim[1]
  y1 <- plot.dim[2]
  
  
  if (plot.type == "confplot") {
    if (missing(main)) main <- paste("Group Configurations") else main <- main
    if (missing(xlab)) xlab <- paste("Dimension", x1,sep = " ") else xlab <- xlab
    if (missing(ylab)) ylab <- paste("Dimension", y1,sep = " ") else ylab <- ylab

    if (missing(xlim)) xlim <- range(x$gspace[,x1])*1.1
    if (missing(ylim)) ylim <- range(x$gspace[,y1])*1.1
  
        
    plot(x$gspace[,x1], x$gspace[,y1], main = main, type = type, xlab = xlab, ylab = ylab, 
         xlim = xlim, ylim = ylim, pch = pch, asp = asp, col = col, ...)
    if (label.conf$label) text(x$gspace[,x1], x$gspace[,y1], labels = rownames(x$gspace), 
                               cex = label.conf$cex, pos = label.conf$pos, 
                               col = label.conf$col)
    
    if (identify) {
       identify(x$gspace[,x1], x$gspace[,y1], labels = rownames(x$gspace), cex = 0.8)
    }
  }

  #----------------------- Stress decomposition -----------------
  if (plot.type == "stressplot") {
    if (missing(main)) main <- paste("Stress Decomposition Chart") else main <- main
    if (missing(xlab)) xlab <- "Objects" else xlab <- xlab
    if (missing(ylab)) ylab <- "Stress Proportion (%)" else ylab <- ylab

    spp.perc <- sort(x$spp, decreasing = TRUE)
    xaxlab <- names(spp.perc)

    if (missing(xlim)) xlim1 <- c(1,length(spp.perc)) else xlim1 <- xlim
    if (missing(ylim)) ylim1 <- range(spp.perc) else ylim1 <- ylim
    
    plot(1:length(spp.perc), spp.perc, xaxt = "n", type = "p",
         xlab = xlab, ylab = ylab, main = main, xlim = xlim1, ylim = ylim1, ...)
    text(1:length(spp.perc), spp.perc, labels = xaxlab, pos = 3, cex = 0.8)
    for (i in 1:length(spp.perc)) lines(c(i,i), c(spp.perc[i],0), col = "lightgray", lty = 2)
                                  
  }

  #------------------------------ bubble plot -------------------------
  if (plot.type == "bubbleplot")
  {

    if (missing(main)) main <- paste("Bubble Plot") else main <- main
    if (missing(xlab)) xlab <- paste("Dimension", x1,sep = " ") else xlab <- xlab
    if (missing(ylab)) ylab <- paste("Dimension", y1,sep = " ") else ylab <- ylab

    if (missing(xlim)) xlim <- range(x$gspace[,x1])*1.1
    if (missing(ylim)) ylim <- range(x$gspace[,y1])*1.1
    
    spp.perc <- x$spp/sum(x$spp)*100
    bubsize <- spp.perc/length(spp.perc)*(bubscale + 3)
    
    plot(x$gspace, cex = bubsize, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim)
    xylabels <- x$gspace
    ysigns <- sign(x$gspace[,y1])
    xylabels[,2] <- (abs(x$gspace[,y1])-(x$gspace[,y1]*(bubsize/50)))*ysigns 
    text(xylabels, rownames(x$gspace), pos = 3,cex = 0.7)    
  }  

}
