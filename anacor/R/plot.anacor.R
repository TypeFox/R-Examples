`plot.anacor` <-
function(x, plot.type = "jointplot", plot.dim = c(1,2), legpos = "top", arrows = FALSE, conf = 0.95, 
         wlines = 0, xlab, ylab, main, type, xlim, ylim, cex.axis2, ...)
{
# x         ... object of class "anacor"
# plot.type ... string with values "jointplot", "rowplot", "colplot", "benzplot", 
#               "regplot", "transplot", "orddiag"
# plot.dim  ... dimensions to be plotted 
# conf ... confidence ellipsoids, either NULL of conf-value. 

options(locatorBell = FALSE)
n <- dim(x$tab)[1]
m <- dim(x$tab)[2]
if (x$ndim < 2) stop("No 2D plot can be produced for ndim < 2!")

if ((plot.type == "benzplot") && (!any(x$scaling == "Benzecri"))) 
  stop("A benzplot can be produced for Benzecri scaled scores only!")

#------------------------------------ regplot -----------------------------------
#draws a before and after regression plot for a table and a set of scores

if (plot.type == "regplot") {
  if (length(plot.dim) != 1) plot.dim <- plot.dim[1]
  if (missing(xlab)) xlab1 = "row" else xlab1 <- xlab
  if (missing(ylab)) ylab1 = "column" else ylab1 <- ylab
  if (missing(main)) main1 <- paste("Unscaled Solution Dimension",plot.dim) else main1 <- main
  if (missing(main)) main2 <- paste("Scaled Solution Dimension",plot.dim) else main2 <- main
  if (missing(type)) type <- "b"
    
  #before ca plot
  tau <- sum(x$tab)
  pr <- x$tab/tau                                     #relative frequencies
  r <- rowSums(pr)                                    #relative row margin
  c <- colSums(pr)                                    #relative column margins
  xave <- as.vector(as.matrix(pr)%*%1:m)/r            #
  yave <- as.vector(1:n%*%as.matrix(pr))/c
  z <- c(1:n,1:m) 
  plot(z, z, type = "n", xlab = paste(xlab1," categories"), ylab = paste(ylab1, " categories"), main = main1, 
  xaxt = "n", yaxt = "n", xlim = c(1,n), ylim = c(1,m), ...)
  axis(1, at = 1:(dim(x$tab)[1]), labels = rownames(x$tab), ...)
  axis(2, at = 1:(dim(x$tab)[2]), labels = colnames(x$tab), ...)
  points(1:n, xave, type = type, col = "RED")
  points(yave, 1:m, type= type, col = "BLUE")
  abline(v=1:n, h=1:m, col = "lightgray", lty = 2 )
  for (i in 1:n) text(rep((1:n)[i],m),1:m,as.character(x$tab[i,]),cex=.8, col = "lightgray")

 
  dev.new()
  #scaled solution
  
  scalevec <- c(as.vector(x$row.scores[,plot.dim]),as.vector(x$col.scores[,plot.dim]))
  if (missing(xlim)) xlim <- range(scalevec)
  if (missing(ylim)) ylim <- range(scalevec)
  
  xa <- x$row.scores[,plot.dim]
  ya <- x$col.scores[,plot.dim]
  xave <- as.vector(as.matrix(pr)%*%ya)/r
  yave <- as.vector(xa%*%as.matrix(pr))/c
  z <- c(xa,ya) 
  plot(z, z, type = "n", xlab = paste(xlab1," scores"), ylab = paste(ylab1," scores"),main = main2,  
  xlim = xlim, ylim = ylim,...)
  
  points(xa, xave, type = type, col = "RED")
  points(yave, ya, type = type, col = "BLUE")
  abline(v = xa, h = ya, col = "lightgray", lty = 2)

  if (missing(cex.axis2)) cex.axis2 <- 0.6
   
  for (i in 1:n) text(rep(xa[i],m),ya,as.character(x$tab[i,]), col = "lightgray")

  #axis(3, at = xa, labels = names(xa), cex.axis = 0.6, col.axis = "lightgray", padj = 1,...)
  #axis(4, at = ya, labels = names(ya), cex.axis = 0.6, col.axis = "lightgray", padj = -1,...)
  
  axis(3, at = xa, labels = names(xa), cex.axis = cex.axis2, col.axis = "lightgray", padj = 1,...)
  axis(4, at = ya, labels = names(ya), cex.axis = cex.axis2, col.axis = "lightgray", padj = -1,...)
  
}

#-------------------------------- end regplot -------------------------------


#-------------------------------- transplot ------------------------------------
#draws row and column or row/column transformation plots.

if (plot.type == "transplot") {
  #if (length(plot.dim) != 1) plot.dim <- plot.dim[1]
  if (missing(main)) main1 <- paste("Transformation Plot - Rows") else main1 <- main
  if (missing(main)) main2 <- paste("Transformation Plot - Columns") else main2 <- main
  if (missing(type)) type <- "b"
  
  xa <- x$row.scores[,plot.dim]
  ya <- x$col.scores[,plot.dim]
    
  matplot(1:n, xa, type = type, xlab = "row categories", ylab = "row scores", main = main1, 
  xaxt = "n", pch = 1,...)
	axis(1, at = 1:(dim(x$tab)[1]), labels = rownames(x$tab))
	abline(v = 1:n, col = "lightgray", lty = 2)
	legend(legpos, legend = paste("Dimension ",plot.dim), lty = 1:length(plot.dim), 
  col = 1:length(plot.dim), bg = "white", cex = 0.8)
  
  dev.new()
  matplot(1:m, ya, type = type, xlab = "column categories", ylab = "column scores", 
  main = main2, xaxt = "n", pch = 1,...)
  axis(1, at = 1:(dim(x$tab)[2]), labels = colnames(x$tab))
  abline(v = 1:m, col = "lightgray", lty = 2)
  legend(legpos, legend = paste("Dimension ",plot.dim), lty = 1:length(plot.dim),
  col = 1:length(plot.dim), bg = "white", cex = 0.8)
}

#--------------------------------- end transplot -------------------------------


#----------------------------------- rowplot -----------------------------------
if (plot.type == "rowplot") {
  if (length(plot.dim) != 2) stop("plot.dim must be of length 2")
  if (missing(main)) main1 <- "Row Plot" else main1 <- main
  if (missing(xlab)) xlab <- paste("Dimension",plot.dim[1])
  if (missing(ylab)) ylab <- paste("Dimension",plot.dim[2])
    
  xa <- x$row.scores[,plot.dim]
  if (missing(xlim)) xlim <- range(xa)*1.2
  if (missing(ylim)) ylim <- range(xa)*1.2

  plot(xa, type = "p", xlab = xlab, ylab = ylab, main = main1, pch = 20, xlim = xlim, ylim = ylim,...)
  text(xa, rownames(x$row.scores), pos = 3, cex = 0.8, offset = 0.2, ...)
  if ((!is.null(conf)) && (!is.null(x$row.acov))) {
    rad<-sqrt(qchisq(conf,2))
    for (i in 1:n) ellipse(xa[i,], x$row.acov[plot.dim,plot.dim,i], rad, col = "darkgray", center.cex = 0.4, 
    lwd = 1)
  }
  if (arrows) arrows(0,0,xa[,1],xa[,2], length = 0.1)
  abline(h = 0, v = 0, col = "lightgray", lty = 2)
 
  
}
#---------------------------------- end rowplot --------------------------------


#----------------------------------- colplot -----------------------------------
if (plot.type == "colplot") {
  if (length(plot.dim) != 2) stop("plot.dim must be of length 2")
  if (missing(main)) main1 <- "Column Plot" else main1 <- main
  if (missing(xlab)) xlab <- paste("Dimension",plot.dim[1])
  if (missing(ylab)) ylab <- paste("Dimension",plot.dim[2])
    
  ya <- x$col.scores[,plot.dim]
  if (missing(xlim)) xlim <- range(ya)*1.2
  if (missing(ylim)) ylim <- range(ya)*1.2
  
  plot(ya, type = "p", xlab = xlab, ylab = ylab, main = main1, pch = 20, xlim = xlim, ylim = ylim,...)
  text(ya, rownames(x$col.scores), pos = 3, cex = 0.8, offset = 0.2)
  if ((!is.null(conf)) && (!is.null(x$col.acov))) {
    rad<-sqrt(qchisq(conf,2))
    for (i in 1:m) ellipse(ya[i,], x$col.acov[plot.dim,plot.dim,i], rad, col = "darkgray", center.cex = 0.4, 
    lwd = 1)
  }
  if (arrows) arrows(0,0,ya[,1],ya[,2], length = 0.1)
  abline(h = 0, v = 0, col = "lightgray", lty = 2)
}
	
#------------------------------- end colplot -----------------------------------

#--------------------------------- jointplot -----------------------------------
#jointplot makes a labelled biplot, using one of four different scalings.
#graph plot is optionally added 

if (plot.type == "jointplot") {
  if (length(plot.dim) != 2) stop("plot.dim must be of length 2")
  if (missing(xlab)) xlab <- paste("Dimension",plot.dim[1])
  if (missing(ylab)) ylab <- paste("Dimension",plot.dim[2])
  if (missing(main)) main1 <- paste("Joint plot") else main1 <- main 
  
  xa <- x$row.scores[,plot.dim]
  ya <- x$col.scores[,plot.dim]

  if (missing(xlim)) xlim <- range(rbind(xa,ya))*1.2
  if (missing(ylim)) ylim <- range(rbind(xa,ya))*1.2                                 
  
  plot(xa, type = "p", xlab = xlab, ylab = ylab, main = main1, pch = 20, xlim = xlim, 
  ylim = ylim, col = "RED",...)
  points(ya, pch = 20, col = "BLUE")
  text(xa, rownames(xa), col = "RED", pos = 3, cex = 0.8, offset = 0.2)
  text(ya, rownames(ya), col = "BLUE", pos = 3, cex = 0.8, offset = 0.2)
  if (!is.null(conf)) {
    rad<-sqrt(qchisq(conf,2))
    if (!is.null(x$row.acov)) 
    {
    for (i in 1:n) ellipse(xa[i,], x$row.acov[plot.dim,plot.dim,i], rad, col = hcl(30), center.cex = 0.4, lwd = 1)
    }
    if (!is.null(x$row.acov)) {
    for (i in 1:m) ellipse(ya[i,], x$col.acov[plot.dim,plot.dim,i], rad, col = hcl(230), center.cex = 0.4, lwd = 1)
    }
  }
  if (arrows) {
    arrows(0,0,xa[,1],xa[,2], col = hcl(30), length = 0.1)
    arrows(0,0,ya[,1],ya[,2], col = hcl(230), length = 0.1)
  }  
  abline(h = 0, v = 0, col = "lightgray", lty = 2)
}

#--------------------------------- end jointplot -------------------------------

#----------------------------------- graphplot ---------------------------------
if (plot.type == "graphplot") {
  if (length(plot.dim) != 2) stop("plot.dim must be of length 2")
  if (missing(xlab)) xlab <- paste("Dimension",plot.dim[1])
  if (missing(ylab)) ylab <- paste("Dimension",plot.dim[2])
  if (missing(main)) main <- paste("Graphplot")  
  
  xa <- x$row.scores[,plot.dim]
  ya <- x$col.scores[,plot.dim]
  if (missing(xlim)) xlim <- range(rbind(xa,ya))
  if (missing(ylim)) ylim <- range(rbind(xa,ya))
  
  
  plot(xa, type = "p", xlab = xlab, ylab = ylab, main = main, pch = 16, col = "RED",
  xlim = xlim, ylim = ylim)
  points(ya, pch = 2, col = "BLUE")
  
  rg <- round(wlines*x$tab/max(x$tab))
  if (wlines > 0) {
  	for (i in 1:n) for (j in 1:m)
  		if (rg[i,j] > 0) lines(c(xa[i,1],ya[j,1]),c(xa[i,2],ya[j,2]),lwd=rg[i,j])
  } else {
    for (i in 1:n) for (j in 1:m)
  		lines(c(xa[i,1],ya[j,1]),c(xa[i,2],ya[j,2]))
  }
  identify(rbind(xa,ya),labels = c(rownames(x$tab),colnames(x$tab)), cex = 0.8)
}

#------------------------------ end graphplot-----------------------------------

#----------------------------------- benzplot ----------------------------------
if (plot.type == "benzplot") 
{
  if(x$scaling[1] == "Benzecri")
  {
    do <- x$bdmat[[1]]                        #observed row distances
    dr <- x$bdmat[[2]]                        #fitted row distances
    dmax <- max(dr,do)
    idlabels.mat <- expand.grid(rownames(do),rownames(dr))
    idlabels <- apply(idlabels.mat, 1, function(xx) {
                                        paste(xx[1], "-", xx[2])
                                       })
    if (missing(main)) main1 <- paste("Benzecri Distances - Rows") else main1 <- paste(main," - Rows")

    plot(do, dr, xlim = c(0,dmax), ylim = c(0,dmax),  xlab = "observed distances", ylab = "fitted distances", 
       main = main1, axes = FALSE, type = "p", col = "RED",...)
    abline(0,1)                               #diagonal
    abline(h = 0, v = 0)                      #axes
    for (i in 1:n)
  	for (j in 1:i)
            lines(c(do[i,j],do[i,j]),c(0,dr[i,j]), col = "lightgray", lty = 2)
    identify(cbind(as.vector(do),as.vector(dr)), labels = idlabels, col = "RED")
  }

  if(x$scaling[2] == "Benzecri")
  {
   do <- x$bdmat[[3]]                          #observed column distances
   dr <- x$bdmat[[4]]                          #observed row distances
   idlabels.mat <- expand.grid(rownames(do),rownames(dr))
   idlabels <- apply(idlabels.mat, 1, function(xx) {
                                        paste(xx[1], "-", xx[2])
                                       })
   dmax <- max(dr,do)
   if (missing(main)) main2 <- paste("Benzecri Distances - Columns") else main2 <- paste(main," - Columns")
  
   dev.new()
   plot(do, dr, xlim = c(0,dmax), ylim = c(0,dmax), xlab = "observed distances", ylab = "fitted distances",
        main = main2, axes=FALSE, col = "BLUE",...)
   abline(0,1)
   abline(h = 0, v = 0)
   for (i in 1:m)
   	for (j in 1:i)
   		lines(c(do[i,j],do[i,j]),c(0,dr[i,j]), col = "lightgray", lty = 2)
   identify(cbind(as.vector(do),as.vector(dr)), labels = idlabels, col = "BLUE")
  } 
}
#----------------------------------end benzplot --------------------------------

#-------------------------------- ordination diagram ---------------------------
if (plot.type == "orddiag") {
  if (length(plot.dim) != 2) stop("plot.dim must be of length 2")
  if (missing(xlab)) xlab <- paste("Dimension",plot.dim[1])
  if (missing(ylab)) ylab <- paste("Dimension",plot.dim[2])
  if (missing(main)) main1 <- paste("Ordination Diagram") else main1 <- main 
  
  xa <- x$row.scores[,plot.dim]
  ya <- x$col.scores[,plot.dim]
 
  crcor.t <- do.call(rbind, x$isetcor)
  if (is.null(crcor.t)) stop("Ordination Diagram works for CCA only!") else crcor <- t(crcor.t)

  if (missing(xlim)) xlim <- range(rbind(xa,ya, crcor))
  if (missing(ylim)) ylim <- range(rbind(xa,ya, crcor))                                 
  
  plot(xa, type = "p", xlab = xlab, ylab = ylab, main = main1, pch = 20, xlim = xlim, 
  ylim = ylim, col = "RED",...)
  points(ya, pch = 20, col = "BLUE")
  text(xa, rownames(xa), col = "RED", pos = 3, cex = 0.8, offset = 0.2)
  text(ya, rownames(ya), col = "BLUE", pos = 3, cex = 0.8, offset = 0.2)
  abline(h = 0, v = 0, col = "lightgray", lty = 2)

  #points(crcor, pch = 20, col = "GRAY")
  text(crcor, rownames(crcor), col = "GRAY", cex = 0.8)
  arrows(0,0,crcor[,1],crcor[,2], length = 0.05, col = "GRAY")

}


#-------------------------------- end ordination diagram ------------------------
}

