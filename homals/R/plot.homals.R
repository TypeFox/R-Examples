`plot.homals` <-
function(x, plot.dim = c(1,2), plot.type = "loadplot", var.subset, main, type, xlab, ylab, 
         xlim, ylim, leg.pos = "topright", identify = FALSE, ...)
{
#S3 plot method for objects of class "homals"
#Produces various 2D-plots
#plot.dim ... vector of length 2 with dimensions to be plotted against
#plot.type ... type of plot to be drawn: "catplot","graphplot","hullplot","labplot",
#              "lossplot","objplot","prjplot","spanplot","starplot", "trfplot", 
#              "vecplot","vorplot", "jointplot","loadplot".
#var.subset ... numeric vector with subset of variables




#plot.type <- plot.type[1]          #use first plot-type only

options(locatorBell = FALSE)
if (x$ndim == 1) stop("No plots can be drawn for ndim = 1 !")
if (length(plot.dim) !=  2) stop("plot.dim must be of length 2!")
if ((plot.type != "trfplot") && (plot.type != "screeplot")) {      #plot.dim are ignored for trfplot
  pd1 <- plot.dim[1]
  pd2 <- plot.dim[2]
  if (pd2 > x$ndim) stop("Only",x$ndim,"dimensions were extracted!")
  if (missing(xlab)) xlab <- paste("Dimension",pd1)
  if (missing(ylab)) ylab <- paste("Dimension",pd2)
}

nvar <- dim(x$dframe)[2]
if (missing(var.subset)) var.subset <- 1:nvar
   

#----------------------------------loadplot-------------------------------------
if (plot.type == "loadplot") {
  xycoor <- t(sapply(x$loadings, function(xy) xy[1,c(pd1,pd2)]))
  if (missing(main)) main1 <- "Loadings plot" else main1 <- main

  xlim.min <- min(xycoor[,1],0)
  xlim.max <- max(xycoor[,1],0)
  ylim.min <- min(xycoor[,2],0)
  ylim.max <- max(xycoor[,2],0)
  if (missing(xlim)) xlim <- c(xlim.min,xlim.max)*1.2
  if (missing(ylim)) ylim <- c(ylim.min,ylim.max)*1.2

   
  plot(xycoor,type = "p", pch = 20, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main1, cex = 0.5)
  for (i in 1:nvar) arrows(0, 0, xycoor[i,1],xycoor[i,2], length = 0.08)   #lines(rbind(xycoor[i,],c(0,0)))
  abline(h = 0, col = "lightgray", lty = 2)
  abline(v = 0, col = "lightgray", lty = 2)
  
  if (identify) {
    identify(xycoor, labels = rownames(xycoor), cex = 0.7)
  } else {
    posvec <- apply(xycoor, 1, sign)[2,] + 2      
    text(xycoor, labels = rownames(xycoor), pos = posvec, cex = 0.7)
  }  
}
#-------------------------------- end loadplot ---------------------------------




#----------------------------------catplot--------------------------------------
#plots the rank-restricted category quantifications for each variable

if (plot.type == "catplot") {

  if (missing(type)) type <- "b"
  if (missing(xlim)) xlim <- range(sapply(x$catscores, function(zz) range(zz[,c(pd1,pd2)])))
  if (missing(ylim)) ylim <- range(sapply(x$catscores, function(zz) range(zz[,c(pd1,pd2)])))

  for (i in var.subset) {
    if (missing(main)) main1 <- paste("Category plot for",colnames(x$dframe[i]))  else main1 <- main

    par("ask"=TRUE)
    plot(x$catscores[[i]][,c(pd1,pd2)], type = type, xlim = xlim, ylim = ylim,
    main = main1, xlab = xlab, ylab = ylab, ...)
    text(x$catscores[[i]][,c(pd1,pd2)], levels(x$dframe[,i]), pos = 3, cex = 0.7)
    abline(h = 0, v = 0, col = "gray", lty = 2)
  }
}
#----------------------------------end catplot----------------------------------

#------------------------------------jointplot----------------------------------
if (plot.type == "jointplot") {
  xylist <- lapply(x$catscores, apply, 2, range)         
  xytab <- sapply(xylist, function(yy) yy[,c(pd1,pd2)])
  xmin <- min(xytab[1,], x$objscores[,pd1])
  xmax <- max(xytab[2,], x$objscores[,pd1])
  ymin <- min(xytab[3,], x$objscores[,pd2])
  ymax <- max(xytab[4,], x$objscores[,pd2])
    
  #xylim <- c(as.vector(xytab), as.vector(x$objscores))
  if (missing(xlim)) xlim <- c(xmin, xmax)
  if (missing(ylim)) ylim <- c(ymin, ymax)
  
  if (missing(main)) main <- "Joint Plot"
  plot(x$objscores[,c(pd1,pd2)], type = "n", main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...) 
  text(x$objscores[,c(pd1,pd2)], labels = rownames(x$dframe), col = 1)
  catcol <- rainbow(ncol(x$dframe))
  
  catleg <- NULL
  for (j in var.subset)
  {
    text(x$catscores[[j]][,c(pd1,pd2)], labels = rownames(x$catscores[[j]]), col = catcol[j], cex = 0.8)
    catleg <- c(catleg, catcol[j])
  }
  legend(leg.pos,colnames(x$dframe)[var.subset], col = catleg, pch = 22)
}

#----------------------------------end jointplot--------------------------------

#------------------------------------graphplot----------------------------------
if (plot.type == "graphplot") {

  if (missing(main)) main <- "Graphplot"
  if (missing(xlim)) xlim <- range(x$objscores[,pd1])*1.2
  if (missing(ylim)) ylim <- range(x$objscores[,pd2])*1.2
  plot(x$objscores[,c(pd1,pd2)], col = "GREEN", pch = 8, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,...)           #draw scores

  dmat <- NULL
  for (j in 1:ncol(x$dframe))
  {
    y <- computeY(x$dframe[,j], x$objscores[,c(pd1,pd2)])
    dmat <- rbind(dmat, y)
    points(y, col = "RED", pch = 16)                             #insert points
    for (i in 1:nrow(x$dframe))
      lines(rbind(x$objscores[i,c(pd1,pd2)], y[x$dframe[i,j],]))                #insert lines
  }
  repvec <- sapply(x$catscores, function(yy) dim(yy)[1])
  varnames <- rep(colnames(x$dframe),repvec)
  rownames(dmat) <- paste(varnames, rownames(dmat))
  xycoor <- rbind(dmat, x$objscores[,c(pd1,pd2)])
  if (identify) { identify(xycoor, labels = rownames(xycoor), cex = 0.7)
  } else { text(xycoor, labels = rownames(xycoor), cex = 0.7, pos = 3) }
}

#----------------------------------end graphplot--------------------------------

#------------------------------------hullplot-----------------------------------
#plots the convex hulls
if (plot.type == "hullplot") {

  if (missing(xlim)) xlim <- range(x$objscores[,pd1])
  if (missing(ylim)) ylim <- range(x$objscores[,pd2])

  for (i in var.subset) {
    
    if (missing(main)) main1 <- paste("Hullplot for",colnames(x$dframe[i]))  else main1 <- main
    
    par("ask" = TRUE) 
    plot(x$objscores[,c(pd1,pd2)], col = "GREEN", pch = 8, main = main1, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
       
    for (j in levels(x$dframe[,i])) 
    {
      ind <- which(j==x$dframe[,i])                      #object index for convex hulls
  	  lst <- ind[chull(x$objscores[ind,c(pd1,pd2)])]                  #convex hull over ind
  	  lines(x$objscores[c(lst,lst[1]),c(pd1,pd2)])
  	  text(x$objscores[lst,c(pd1,pd2)], j, cex = 0.7)
    }
  }
}
#-----------------------------------end hullplot--------------------------------

#--------------------------------------labplot----------------------------------
#plot labeled object scores (for each variable separately)

if (plot.type == "labplot") {

  if (missing(xlim)) xlim <- range(x$objscores[,pd1])
  if (missing(ylim)) ylim <- range(x$objscores[,pd2])
  
  for (i in var.subset) {
    if (missing(main)) main1 <- paste("Labplot for",colnames(x$dframe[i]))  else main1 <- main
    par("ask" = TRUE)
    plot(x$objscores[,c(pd1,pd2)], type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main1, ...)
    text(x$objscores[,c(pd1,pd2)], as.vector(x$dframe[,i]), cex = 0.7)
  }
} 
#-----------------------------------end labplot---------------------------------

#------------------------------------lossplot-----------------------------------
if (plot.type == "lossplot") {
  
  for (i in var.subset) { 
    
    if (missing(main)) main1 <- paste("Lossplot for",colnames(x$dframe[i])) else main1 <- main
    
    z <- computeY(x$dframe[,i], x$objscores[,c(pd1,pd2)])
    k <- dim(z)[1]
    
    if (missing(xlim)) xlim1 <- range(c(z[,2],x$catscores[[i]][,pd1])) else xlim1 <- xlim
    if (missing(ylim)) ylim1 <- range(c(z[,2],x$catscores[[i]][,pd2])) else ylim1 <- ylim
    
    par("ask" = TRUE)
    plot(x$catscores[[i]][,c(pd1,pd2)], type = "p", main = main1, xlab = xlab, ylab = ylab,
         xlim = xlim1, ylim = ylim1, col = "RED", ...)
    text(x$catscores[[i]][,c(pd1,pd2)], levels(x$dframe[,i]), pos = 3, col = "RED", cex = 0.7)
    lines(x$catscores[[i]][,c(pd1,pd2)], col = "RED")
    
    points(z,type="p", col = "blue")
    text(z, levels(x$dframe[,i]), col="blue", pos = 3, cex = 0.7)
    lines(z, col="blue")
    for (j in 1:k) lines(rbind(x$catscores[[i]][j,c(pd1,pd2)],z[j,]),col="lightgray", lty=3)
    
    abline(h = 0, v = 0, col = "gray", lty = 2)
  }  
}
#----------------------------------end lossplot---------------------------------

#------------------------------------objplot------------------------------------
#draws labeled object score plot

if (plot.type == "objplot") {
    
  if (missing(xlim)) xlim <- range(x$objscores[,pd1])
  if (missing(ylim)) ylim <- range(x$objscores[,pd2])
  if (missing(main)) main1 <- "Plot Object Scores" else main1 <- main

  plot(x$objscores[,c(pd1,pd2)], type = "n", main = main1, xlab = xlab, ylab = ylab, 
       xlim = xlim, ylim = ylim, ...)
  text(x$objscores[,c(pd1,pd2)], rownames(x$objscores), cex = 0.7)
} 

#---------------------------------end objplot-----------------------------------

#----------------------------------prjplots-------------------------------------
#draws projection plot

if (plot.type == "prjplot") {

  if (missing(xlim)) xlim <- range(x$objscores[,pd1])
  if (missing(ylim)) ylim <- range(x$objscores[,pd2])
  xylim <- c(min(xlim[1],ylim[1]),max(xlim[2],ylim[2]))
  
  for (i in var.subset) {
  
    if (missing(main)) main1 <- paste("Projection Plot for", colnames(x$dframe[i])) else main1 <- main
     
    a <- x$loadings[[i]][c(pd1,pd2)]
    
    if (x$rank.vec[i] == 1) {
      par("ask" = TRUE)
      plot(x$objscores[,c(pd1,pd2)], type = "p", main = main1, xlab = xlab, ylab = ylab,
           xlim = xylim, ylim = xylim, ...)
      text(x$objscores[,c(pd1,pd2)], as.vector(x$dframe[,i]), pos = 3) 
      text(x$catscores[[i]], levels(x$dframe[,i]), col = "RED", pos = 3)
      slope = a[2]/a[1]
      abline(coef = c(0,slope))
      slope = -a[1]/a[2]
      
      for (j in 1:(dim(x$catscores[[i]])[1]))
  	    abline(coef = c(x$catscores[[i]][j,pd2] - slope*x$catscores[[i]][j,pd1],slope), col = "RED")
  	  for (k in 1:(dim(x$objscores)[1])) {
  	    j <- x$dframe[,i][k]
        icpt <- x$catscores[[i]][j,pd2] - slope*x$catscores[[i]][j,pd1] 
  	    u <-(x$objscores[k,pd1] + slope*(x$objscores[k,pd2]-icpt))/(1+slope^2)  
        lines(rbind(x$objscores[k,c(pd1,pd2)],c(u,icpt+slope*u)), col = "BLUE")
  	  }
      abline(h = 0, v = 0, col = "gray", lty = 2, ...)   
    } else {
      warning("No projection plot for ", colnames(x$dframe[i]),"(rank != 1).", call. = FALSE) 
    }
  }
}
#---------------------------------end prjplot-----------------------------------

#-------------------------------------spanplot----------------------------------
if (plot.type == "spanplot") {

  if (missing(xlim)) xlim <- range(x$objscores[,pd1])
  if (missing(ylim)) ylim <- range(x$objscores[,pd2])
  
  for (i in var.subset) {
    if (missing(main)) main1 <- paste("Span plot for", colnames(x$dframe[i])) else main1 <- main
    
    par("ask" = TRUE)
    plot(x$objscores[,c(pd1,pd2)], col = 1, pch = 21, main = main1, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
    lev <- levels(x$dframe[,i])
    rb<-rainbow(length(lev))
    for (k in lev) {
  	  ind <- which(k==x$dframe[,i])
  	  n <- length(ind)
  	  mm <- mst(dist(x$objscores[ind,c(pd1,pd2)]))
  	  for (j in 1:n) {
  		  jnd <- which(1 == as.vector(mm[j,]))
  		  sapply(jnd, function(r) lines(rbind(x$objscores[ind[j],c(pd1,pd2)], x$objscores[ind[r],c(pd1,pd2)]),
        col = rb[which(lev==k)]))
  		}
  	 legend(leg.pos,paste("Category",lev), col = rb, lty = 1, cex = 0.7)
    } 
  }
}
#----------------------------------end spanplot---------------------------------


#------------------------------------starplot-----------------------------------
if (plot.type == "starplot") {

  if (missing(xlim)) xlim <- range(x$objscores[,pd1])
  if (missing(ylim)) ylim <- range(x$objscores[,pd2])
  
  for (i in var.subset) {
    if (missing(main)) main1 <- paste("Star plot for", colnames(x$dframe[i])) else main1 <- main
    
    par("ask" = TRUE)
    plot(x$objscores[,c(pd1,pd2)], col = "BLUE", main = main1, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
    
    z <- computeY(x$dframe[,i], x$objscores[,c(pd1,pd2)])
    points(z, type="o", pch = 24, col = "RED")
    text(z, levels(x$dframe[,i]), col = "RED", pos = 3)
    for (j in 1:length(x$dframe[,i])) 
      lines(rbind(x$objscores[j,c(pd1,pd2)],z[x$dframe[,i][j],]), col = "BLUE")
    if (identify) { identify(x$objscores[,c(pd1,pd2)], labels = rownames(x$dframe), col = "BLUE", cex = 0.7) 
    } else { text(x$objscores[,c(pd1,pd2)], labels = rownames(x$dframe), col = "BLUE", pos = 3, cex = 0.7) }
  }
}
#----------------------------------end starplot---------------------------------

#------------------------------------vecplot------------------------------------
#draws vector plots

if (plot.type == "vecplot") {
  
  if (missing(xlim)) xlim <- range(x$objscores[,pd1])
  if (missing(ylim)) ylim <- range(x$objscores[,pd2])
  xylim <- c(min(xlim[1],ylim[1]),max(xlim[2],ylim[2]))

  for (i in var.subset) {
    if (missing(main)) main1 <- paste("Vector plot for", colnames(x$dframe[i])) else main1 <- main
     
    if (x$rank.vec[i] == 1) {
      a <- x$loadings[[i]][,c(pd1,pd2)]
      par("ask" = TRUE)
      plot(x$objscores[,c(pd1,pd2)], type = "p", main = main1, xlab = xlab, ylab = ylab,
           xlim = xylim, ylim = xylim, ...)
      text(x$objscores[,c(pd1,pd2)], as.vector(x$dframe[,i]), pos = 3, cex = 0.7) 
      text(x$catscores[[i]][,c(pd1,pd2)], levels(x$dframe[,i]), col = "RED", cex = 0.7)
      slope = a[2]/a[1]
      abline(coef = c(0,slope))
      
      for (j in 1:length(x$dframe[,i])) {
        xs <- xe <- x$objscores[j,c(pd1,pd2)]  
        xe[1] <- (xs[1]+(xs[2]*slope))/(1+(slope^2))
        xe[2] <- slope*xe[1]     
  	    lines(rbind(xs,xe), col = "BLUE")
  	  }
      abline(h = 0, v = 0, col = "gray", lty = 2) 
    } else {
      warning("No vector plot for ", colnames(x$dframe[i]),"(rank != 1).", call. = FALSE) 
    }
  }
}
#----------------------------------end vecplot----------------------------------


#------------------------------------trfplot------------------------------------
#draws transformation plots

if (plot.type == "trfplot") {

   if (missing(type)) type <- "l"
   if (missing(xlab)) xlab <- "original scale"
   if (missing(ylab)) ylab <- "transformed scale"

   for (i in var.subset) {
     
     if (missing(main)) main1 <- paste("Transformation plot for", colnames(x$dframe[i])) else main1 <- main
     if (missing(ylim)) ylim1 <- range(x$low.rank[[i]]) else ylim1 <- ylim
     
     p <- dim(x$low.rank[[i]])[2]                           #number of dimensions
     vlev <- rownames(x$low.rank[[1]])
       
     par("ask" = TRUE)                     #first dimensions
     matplot(x$low.rank[[i]], type = type, main = main1, ylim = ylim1, xlab = xlab, 
            ylab = ylab, xaxt = "n", pch = 20, ...)
     if (p != 1) legend(leg.pos,paste("Solution",1:p),col = 1:p, lty = 1:p)
     axis(1, at = 1:dim(x$low.rank[[i]])[1], labels = rownames(x$low.rank[[i]]))
     
   }
}
#----------------------------------end trfplot----------------------------------

#------------------------------------vorplot------------------------------------
#draws voronoi regions

if (plot.type == "vorplot") {

  for (i in var.subset) {
     
     if (missing(main)) main1 <- paste("Voronoi plot for", colnames(x$dframe[i])) else main1 <- main
     xlim.i <- c(x$objscores[,pd1],x$catscores[[i]][,pd1])
     ylim.i <- c(x$objscores[,pd2],x$catscores[[i]][,pd2])
     if (missing(xlim)) xlim1 <- range(xlim.i) else xlim1 <- xlim
     if (missing(ylim)) ylim1 <- range(ylim.i) else ylim1 <- ylim
   
     par("ask" = TRUE)
     plot(x$objscores[,c(pd1,pd2)], type = "n", main = main1, xlab = xlab, ylab = ylab, 
          xlim = xlim1, ylim = ylim1, ...)
     drawEdges(x$catscores[[i]][,c(pd1,pd2)], far = 1000)
     text(x$objscores[,c(pd1,pd2)], as.vector(x$dframe[,i]))
  }
}
#----------------------------------end vorplot----------------------------------

#---------------------------------- screeplot ----------------------------------
if (plot.type == "screeplot") {

   if (missing(main)) main <- "Scree plot"
   if (missing(xlab)) xlab <- "Dimension"
   if (missing(ylab)) ylab <- "Eigenvalue"

  plot(1:x$ndim, x$eigenvalues, type = "b", xlab = xlab, ylab = ylab, main = main, xaxt = "n", ...)
  axis(1, at = 1:x$ndim, labels = 1:x$ndim)
}  
#-------------------------------- end screeplot --------------------------------
 
#----------------------------------- dmplot ------------------------------------
if (plot.type == "dmplot") {

   if (missing(main)) main <- "Discrminination Measures"
   xylim <- range(x$discrim[,c(pd1, pd2)])
   if (missing(xlim)) xlim <- xylim
   if (missing(ylim)) ylim <- xylim

   plot(x$discrim[,c(pd1,pd2)], type = "p", main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
   for (i in 1:nvar) lines(rbind(x$discrim[i,c(pd1,pd2)],c(0,0)))
   text(x$discrim[,c(pd1,pd2)], rownames(x$discrim), pos = 3)
}  

#-------------------------------- end dmplot -----------------------------------

}

