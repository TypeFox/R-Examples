# agram.R
# Author: Ben Murrell, Dan Murrell and Hugh Murrell
# Adapted from corrgram.r - author: Kevin Wright, 
# Copyright 2006 Kevin Wright, License: GPL2

# The corrgram function was derived from the 'pairs' function.
# Code for plotting ellipses was derived from the ellipse package.
# Additional ideas from the plot.corr function in the 'arm' package.

agram <- function (dataSet, method="pearson", ranking=FALSE, order=FALSE, ...) {
   
  textPanel <- function(x = 0.5, y = 0.5, txt, cex, font, srt) {
    text(x, y, txt, cex=cex, font=font, srt=srt)
  }
  
#   localRank <- function(v) {
#     if (ranking) 
#       return(rwt(v))
#     else
#       return(v)
#   }
#  
#  replaced by ........... if (ranking) dataSet <- rwt(dataSet)
  
  localAxis <- function(side, x, y, xpd, bg, col=NULL, main, oma, ...) {
    ## Explicitly ignore any color argument passed in as
    ## it was most likely meant for the data points and
    ## not for the axis.
    if(side %%2 == 1) Axis(x, side=side, xpd=NA, ...)
    else Axis(y, side=side, xpd=NA, ...)
  }
  
  # Don't pass some arguments on to the panel functions
  localPlot <- function(..., main, oma, font.main, cex.main)
    plot(...)
  localLowerPanel <- function(..., main, oma, font.main, cex.main)
    lower.panel(...)
  localUpperPanel <- function(..., main, oma, font.main, cex.main)
    upper.panel(...)
  
  localDiagPanel <- function(..., main, oma, font.main, cex.main)
    diag.panel(...)
  
  col.corrgram <- function(ncol){
    # Colors to use for the corrgram
    # Red > White > Blue
    colorRampPalette(c("navy","red"))(ncol)
  }
  
  panel.pts <- function(x, y, amat=NULL, cmat=NULL, ...){
    # choose a color
    unexp <- 1 - cmat
    if(unexp > 0.01) 
      nonlin <- max(0,(amat-cmat)/unexp)
    else
      nonlin <- 0.01
    ncol <- 14
    pal <- col.corrgram(ncol)
    col.ind <- as.numeric(cut(nonlin, breaks=seq(from=-0, to=1, length=ncol+1),
                              include.lowest=TRUE))
    col.pts <- pal[col.ind]
    
    # plot.xy(xy.coords(localRank(x),localRank(y)), type="p", col=col.pts,...)
    plot.xy(xy.coords(x,y), type="p", col=col.pts,...)
    box(col="lightgray")
  }
  
  panel.pie <- function(x, y, amat=NULL, cmat=NULL, ...){
    # print(paste(toString(amat),toString(cmat)))
    # Coordinates of box
    usr <- par()$usr
    minx <- usr[1] #min(x, na.rm=TRUE)
    maxx <- usr[2] #max(x, na.rm=TRUE)
    miny <- usr[3] #min(y, na.rm=TRUE)
    maxy <- usr[4] #max(y, na.rm=TRUE)
    # Multiply the radius by .97 so the circles do not overlap
    rx <- (maxx-minx)/2 * .97
    ry <- (maxy-miny)/2 * .97
    centerx <- (minx+maxx)/2
    centery <- (miny+maxy)/2
    
    segments <- 60
    angles <- seq(0,2*pi,length=segments)
    circ <- cbind(centerx + cos(angles)*rx, centery + sin(angles)*ry)
    lines(circ[,1], circ[,2], col='gray30',...)
    
    # use amat to draw the red segment
    segments <- round(60*abs(amat),0) # Watch out for the case with 0 segments
    if(segments>0){
      angles <- seq(pi/2, pi/2+(2*pi* - amat), length=segments)
      circ <- cbind(centerx + cos(angles)*rx, centery + sin(angles)*ry)
      circ <- rbind(circ, c(centerx, centery), circ[1,])
      polygon(circ[,1], circ[,2], col="red")
    }
  
    # now use camt to draw a blue segment   
    segments <- round(60*abs(cmat),0) # Watch out for the case with 0 segments
    if(segments>0){
      angles <- seq(pi/2, pi/2+(2*pi* - cmat), length=segments)
      circ <- cbind(centerx + cos(angles)*rx, centery + sin(angles)*ry)
      circ <- rbind(circ, c(centerx, centery), circ[1,])
      polygon(circ[,1], circ[,2], col="navy")
    }
    
  }
  
  # end of helper functions main code starts now
  
  # rank up front if asked to
  if (ranking) dataSet <- data.frame(rwt(dataSet))
  
  # remember par settings that could be changed
  old.par <- par(no.readonly = TRUE) 
  on.exit(par(old.par))
  
  x <- dataSet
   
  # internal parameters
  # panel = panel.shade
  # lower.panel = panel 
  lower.panel=panel.pts
  # upper.panel = panel
  upper.panel=panel.pie
  diag.panel = textPanel 
  text.panel = textPanel
  label.pos = 0.5 
  label.srt=0
  cex.labels = NULL 
  font.labels = 1
  row1attop = TRUE 
  dir="left" 
  gap = 0

  
  if(is.null(order)) order <- FALSE
    
  if (ncol(x) < 2) stop("need two or more variables (cols) in the input data set")
    
  # Remove non-numeric columns from data frames
  x <- x[ , sapply(x, is.numeric)]

  # x <- rwt(x)
  
  # calculate the association matrix
  amat <- tap(x) 
  
  # make sure amat is positive
  amat <- abs(amat)
  
  # calculate the correlation matrix
  cmat <- cor(x,use="complete.obs",method=method)^2
  
  # Re-order the data to group highly associated variables
  if(order==TRUE | order=="PC" | order=="PCA"){
    # Order by angle size between PCAs (first two) of correlation matrix
    x.eigen <- eigen(amat)$vectors[,1:2]
    e1 <- x.eigen[,1]
    e2 <- x.eigen[,2]
    alpha <- ifelse(e1>0, atan(e2/e1), atan(e2/e1)+pi)
    ord <- order(alpha)
    x <- x[,ord]
    amat <- amat[ord,ord]
    cmat <- cmat[ord,ord]
  } else if (order=="OLO") {
    distx <- dist(amat)
    ss <- seriate(distx, method="OLO")
    ord <- get_order(ss)
    x <- x[,ord]
    amat <- amat[ord,ord]
    cmat <- cmat[ord,ord]
  } else if(order!=FALSE){
    stop("Unknown order argument in 'corrgram'.")
  }
  
  
  
  dots <- list(...); nmdots <- names(dots)
  
 # Get panel functions 
  lower.panel <- match.fun(lower.panel)
  upper.panel <- match.fun(upper.panel)
  diag.panel <- match.fun( diag.panel)
  

  # Plot layout  
  nc <- ncol(x)
  has.labs <- TRUE
  labels <- colnames(x)
  if (is.null(labels)) labels <- paste("var", 1:nc)
  
  oma <- if("oma" %in% nmdots) dots$oma else NULL
  main <- if("main" %in% nmdots) dots$main else NULL
  if (is.null(oma)) {
    oma <- c(4, 4, 4, 4)
    if (!is.null(main)) oma[3] <- 6 # Space for the title
  }
  opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
  on.exit(par(opar))
  
  # Main loop to draw each panel
  for (i in 1:nc)
    for (j in 1:nc) {
      # start a plot but do nothing
      # localPlot(localRank(x[, j]), localRank(x[, i]), xlab = "", ylab = "",
      #           axes = FALSE, type = "n", ...)
      localPlot(x[, j], x[, i], xlab = "", ylab = "",
                axes = FALSE, type = "n", ...)
      # print(paste("i=",toString(i),", j=",toString(j)))
        if(i == j) {
          # Diagonal panel         
          # Diagonal text
          if (has.labs) {
            par(usr = c(0, 1, 0, 1))
            if(is.null(cex.labels)) {
              l.wid <- strwidth(labels, "user")
              cex.labels <- max(0.8, min(2, .9 / max(l.wid)))
            }
            text.panel(0.5, label.pos, labels[i],
                       cex = cex.labels, font = font.labels, srt=label.srt)
          }
        } else if(i < j) {
          # Upper panel
            localUpperPanel(as.vector(x[, j]), as.vector(x[, i]), 
                            amat[j,i], cmat[j,i], ...)
        } else {
          # Lower panel
            localLowerPanel(as.vector(x[, j]), as.vector(x[, i]), 
                            amat[j,i], cmat[j,i], ...)          
        }
      
    }
  
  
  if (!is.null(main)) {
    font.main <- if("font.main" %in% nmdots) dots$font.main else par("font.main")
    cex.main <- if("cex.main" %in% nmdots) dots$cex.main else par("cex.main")
    mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
  }
  

  
  invisible(NULL)
}



# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


