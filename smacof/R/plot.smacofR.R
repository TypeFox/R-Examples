# plot method for rectangular smacof

plot.smacofR <- function(x, plot.type = "confplot", joint = TRUE, plot.dim = c(1,2), 
                         col.rows = hcl(0), 
                         col.columns = hcl(240),
                         label.conf.rows = list(label = TRUE, pos = 3, col = hcl(0, l = 50), cex = 0.8), 
                         label.conf.columns = list(label = TRUE, pos = 3, col = hcl(240, l = 50), cex = 0.8), 
                         type = "p", pch = 20, asp = 1, main, xlab, ylab, xlim, ylim, ...)

# x ... object of class smacofR
# plot.type ... types available: "confplot", "Shepard", "stressplot", "resplot"
# joint ... if TRUE, row and column configurations in 1 plot
  
{
  x1 <- plot.dim[1]
  y1 <- plot.dim[2]
  
  ## --- check label lists
  if (is.null(label.conf.rows$label)) label.conf.rows$label <- TRUE
  if (is.null(label.conf.rows$pos)) label.conf.rows$pos <- 3
  if (is.null(label.conf.rows$col)) label.conf.rows$col <- hcl(0, l = 50)
  if (is.null(label.conf.rows$cex)) label.conf.rows$cex <- 0.8
  if (is.null(label.conf.columns$label)) label.conf.columns$label <- TRUE
  if (is.null(label.conf.columns$pos)) label.conf.columns$pos <- 3
  if (is.null(label.conf.columns$col)) label.conf.columns$col <- hcl(240, l = 50)
  if (is.null(label.conf.columns$cex)) label.conf.columns$cex <- 0.8
  
  
  #--------------------------- configuration plot -----------------------
  if (plot.type == "confplot") {
    #if (missing(type)) type <- "n" else type <- type
    if (missing(xlab)) xlab1 <- paste("Dimension", x1,sep = " ") else xlab1 <- xlab  
    if (missing(ylab)) ylab1 <- paste("Dimension", y1,sep = " ") else ylab1 <- ylab
    ppos.rows <- label.conf.rows[[2]]
    ppos.columns <- label.conf.columns[[2]]
    if (type == "n") ppos.rows <- ppos.columns <- NULL

    if (!joint) {                                                                                 
      op <- par(mfrow = c(1,2))
      if (missing(main)) main1 <- paste("Configuration Plot - Columns") else main1 <- main        #plot column configurations
      #if (missing(xlab)) xlab1 <- paste("Column Configurations D", x1,sep = "") else xlab1 <- xlab
      #if (missing(xlab)) ylab1 <- paste("Column Configurations D", y1,sep = "") else ylab1 <- ylab

      if (missing(xlim)) xlim <- range(x$conf.col[,c(x1,y1)])
      if (missing(ylim)) ylim <- range(x$conf.row[,c(x1,y1)])

      plot(x$conf.col[,x1], x$conf.col[,y1], main = main1, xlab = xlab1, ylab = ylab1, type = type, 
           xlim = xlim, ylim = ylim, col = col.columns, pch = pch, asp = asp, ...)
      if (label.conf.columns[[1]]) text(x$conf.col[,x1], x$conf.col[,y1], labels = rownames(x$conf.col), 
                                        cex = label.conf.columns$cex, pos = ppos.columns, 
                                        col = label.conf.columns[[3]])

      if (missing(main)) main1 <- paste("Configuration Plot - Rows") else main1 <- main       #plot row configurations
      #if (missing(xlab)) xlab1 <- paste("Column Configurations D", x1,sep = "") else xlab1 <- xlab
      #if (missing(ylab)) ylab1 <- paste("Column Configurations D", y1,sep = "") else ylab1 <- ylab
      
      plot(x$conf.row[,x1], x$conf.row[,y1], main = main1, xlab = xlab1, ylab = ylab1, type = type, 
           xlim = xlim, ylim = ylim, col = col.rows, pch = pch, asp = asp, ...)
      if (label.conf.rows[[1]]) text(x$conf.row[,x1], x$conf.row[,y1], labels = rownames(x$conf.row), 
                                     cex = label.conf.rows$cex, pos = ppos.rows, col = label.conf.rows[[3]])      
      par(op)
    
    } else { #joint plot
      if (missing(main)) main1 <- paste("Joint Configuration Plot") else main1 <- main        
      
      fullconf <- rbind(x$conf.col[,c(x1,y1)],x$conf.row[,c(x1,y1)])
      if (missing(xlim)) xlim <- range(fullconf)
      if (missing(ylim)) ylim <- range(fullconf)
      
      plot(x$conf.row[,x1], x$conf.row[,y1], main = main1, xlab = xlab1, ylab = ylab1, type = type, 
             xlim = xlim, ylim = ylim, col = col.rows, pch = pch, asp = asp,...)
      if (label.conf.rows[[1]]) text(x$conf.row[,x1], x$conf.row[,y1], labels = rownames(x$conf.row),
                                     cex = label.conf.rows$cex, pos = ppos.rows, col = label.conf.rows[[3]])    
      points(x$conf.col[,x1], x$conf.col[,y1], main = main1, xlab = xlab1, ylab = ylab1, type = type, 
           xlim = xlim, ylim = ylim, col = col.columns, pch = pch, ...)
      if (label.conf.columns[[1]]) text(x$conf.col[,x1], x$conf.col[,y1], labels = rownames(x$conf.col), 
                                        cex = label.conf.columns$cex, pos = ppos.columns, 
                                        col = label.conf.columns[[3]])
    }
  }

  #---------------- Shepard diagram ------------------
#   if (plot.type == "Shepard") {
#     if (missing(main)) main <- paste("Shepard Diagram") else main <- main
#     if (missing(xlab)) xlab <- "Dissimilarities" else xlab <- xlab
#     if (missing(ylab)) ylab <- "Configuration Distances" else ylab <- ylab
# 
#     if (missing(xlim)) xlim <- range(as.vector(x$obsdiss))
#     if (missing(ylim)) ylim <- range(as.vector(x$confdiss))
# 
#     plot(as.vector(x$obsdiss), as.vector(x$confdiss), main = main, type = "p", pch = 1,
#          xlab = xlab, ylab = ylab, col = "darkgray", xlim = xlim, ylim = ylim, ...)
# 
#     if (!x$metric) {
#       isofit <- isoreg(as.vector(x$obsdiss), as.vector(x$confdiss))  #isotonic regression
#       points(sort(isofit$x), isofit$yf, type = "b", pch = 16)
#     } else {
#       abline(0,1, lty = 2)
#     }
#    
#   }


  #--------------- Residual plot --------------------
#   if (plot.type == "resplot") {
#     if (missing(main)) main <- paste("Residual plot") else main <- main
#     if (missing(xlab)) xlab <- "Configuration Distances" else xlab <- xlab
#     if (missing(ylab)) ylab <- "Residuals" else ylab <- ylab
#     resmat <- residuals(x)
# 
#     fullres <- c(as.vector(x$confdiss), as.vector(resmat)) 
#     if (missing(xlim)) xlim <- range(fullres)
#     if (missing(ylim)) ylim <- range(fullres)
#     
#     plot(as.vector(x$confdiss), as.vector(resmat), main = main, type = "p",
#          xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
#     abline(h = 0, col = "darkgray", lty = 2)  
#   }

  #----------------------- Stress decomposition -----------------
  if (plot.type == "stressplot") {
    if (missing(main)) main1 <- paste("Stress Decomposition Chart - Rows") else main1 <- main
    if (missing(main)) main2 <- paste("Stress Decomposition Chart - Columns") else main2 <- main
    
    if (missing(xlab)) xlab1 <- "Row Objects" else xlab1 <- xlab
    if (missing(xlab)) xlab2 <- "Column Objects" else xlab2 <- xlab
    if (missing(ylab)) ylab <- "Stress Proportion (%)" else ylab <- ylab

   
    #row-wise
    spp.perc.row <- sort((x$spp.row/sum(x$spp.row)*100), decreasing = TRUE)
    xaxlab <- names(spp.perc.row)
   
    if (missing(xlim)) xlim1 <- c(1,length(spp.perc.row)) else xlim1 <- xlim
    if (missing(ylim)) ylim1 <- range(spp.perc.row) else ylim1 <- ylim
    
    op <- par(mfrow = c(1,2))
    plot(1:length(spp.perc.row), spp.perc.row, xaxt = "n", type = "p",
         xlab = xlab1, ylab = ylab, main = main1, xlim = xlim1, ylim = ylim1, ...)
    text(1:length(spp.perc.row), spp.perc.row, labels = xaxlab, pos = 3, cex = 0.8)
    for (i in 1:length(spp.perc.row)) lines(c(i,i), c(spp.perc.row[i],0), col = "lightgray", lty = 2)                               
       
    spp.perc.col <- sort((x$spp.col/sum(x$spp.col)*100), decreasing = TRUE)
    xaxlab <- names(spp.perc.col)
   
    if (missing(xlim)) xlim1 <- c(1,length(spp.perc.col)) else xlim1 <- xlim
    if (missing(ylim)) ylim1 <- range(spp.perc.col) else ylim1 <- ylim
    
    plot(1:length(spp.perc.col), spp.perc.col, xaxt = "n", type = "p",
         xlab = xlab2, ylab = ylab, main = main2, xlim = xlim1, ylim = ylim1, ...)
    text(1:length(spp.perc.col), spp.perc.col, labels = xaxlab, pos = 3, cex = 0.8)
    for (i in 1:length(spp.perc.col)) lines(c(i,i), c(spp.perc.col[i],0), col = "lightgray", lty = 2)  
    par(op)  
  }
}
