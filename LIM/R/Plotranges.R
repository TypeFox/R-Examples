
Plotranges <- function(...) UseMethod ("Plotranges")

##==============================================================================
## Plots ranges and names                               ##
##==============================================================================

Plotranges.double <- function (min, max,  value = NULL,
    labels=NULL, log="", pch = 16, pch.col="black",
    line.col = "gray", seg.col="black",
    xlim = NULL, main = NULL, xlab = NULL, ylab = NULL,
    lab.cex = 1.0, mark = NULL, ...)  {


  ##-----------------------------------------------------------------
  ## constructing the data
  ##-----------------------------------------------------------------
  if (! is.vector(value))
    value<-as.vector(unlist(value))
  ranges   <- cbind(min,max,value)
  if (log =="x") {

	  minflow<-min(ranges[ranges!=0])     ## minimum, different from 0
    ranges[ranges==0] <- minflow
    min[min==0]       <- minflow        ## replace 0 with minimum
    max[max==0]       <- minflow
    value[value==0]   <- minflow
  }
  numflows <- length(min)

  ##-----------------------------------------------------------------
  if (is.null(labels)) labels <- names(min)
  if (is.null(labels)) labels <- names(max)
  if (is.null(labels)) labels <- as.character(1:numflows)

  labelwidth   <- max(strwidth(labels, "inch")*lab.cex, na.rm = TRUE)
  labelheight  <- strheight("M", "inch")
  plot.new()

  ##-----------------------------------------------------------------
  ## new margins
  ##-----------------------------------------------------------------

  nmar         <- nm <- par("mar")
  nmar[2]      <- nmar[4] + (labelwidth + 0.1)/labelheight
  par(mar = nmar)

  y            <- 1:numflows
  if (is.null (xlim))
    xlim <- range(ranges, na.rm=TRUE)
  ylim <- c(0, numflows + 1)

  plot.window(xlim = xlim, ylim = ylim, log = log)

  loffset <- (labelwidth + 0.1)/labelheight
  labs    <- labels[y]
  mtext(labs, side = 2, line = loffset, at = y, adj = 0,
        las = 2, cex = par("cex") * lab.cex, ...)

  abline  (h = y, lty = "dotted", col = line.col)
  if(!is.null(value))
    points  (value, y, pch = pch, col = pch.col)
  segments(min,y,max,y,col=seg.col,lty=1)
  if (! is.null(mark)) {
    text(labels=rep("*",length(mark)),
         rep(xlim[2],length(mark)),mark)
  }
  axis(1)
  box()
  title(main = main, xlab = xlab, ylab = ylab, ...)
  invisible()
  par("mar"=nm)

}


##==============================================================================
## plots ranges and names for LIM input
##==============================================================================

Plotranges.lim <- function (lim=NULL, labels=NULL, type="X", log="",
   pch = 16, pch.col="black", line.col = "gray", seg.col="black",
   xlim = NULL, main = NULL, xlab = NULL, ylab = NULL,
   lab.cex = 1.0, index=NULL, ...)   {

  ##-----------------------------------------------------------------
  ## constructing the data
  ##-----------------------------------------------------------------
  if (type == "V") {
    ranges <- Varranges(lim)
    value  <- Variables(lim)

  } else {
    ranges <- Xranges(lim)
    value  <- Lsei.lim (lim, parsimonious=TRUE)$X
  }
  value <- unlist(value)
  if (is.null(index)) {
    index <- 1:nrow(ranges)
  } else {
    ranges <- ranges[index,]
    value <- value[index]
  }

  infinity <- which (ranges[,2]==1e30)
  if (length(infinity) >0)
    ranges[infinity,2]<- NA else infinity <- NULL

  if (is.null(labels))
    labels <- rownames(ranges)
  if (is.null(labels))
    labels <- as.character(1:length(value))

  Plotranges.double(min=ranges[,1],max=ranges[,2],value =value,
     labels=labels,log=log,pch=pch,pch.col=pch.col,line.col=line.col,
     seg.col=seg.col,xlim=xlim,main=main,xlab=xlab,ylab=ylab,
     lab.cex=lab.cex,mark=infinity,...)

}

##==============================================================================
## plots ranges and names for file input
##==============================================================================

Plotranges.character <- function (file, ...)  {

  ##-----------------------------------------------------------------
  ## constructing the data
  ##-----------------------------------------------------------------

  lim <- Setup(file)
  Plotranges.lim(lim,...)

}

