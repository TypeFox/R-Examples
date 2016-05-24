eccomp <- function(ft, xlim, ylim, main, xlab, ylab, do.points=TRUE,
                    datapch, datacol, fitlty, fitcol, addlegend = TRUE, 
                   legendtext, xlegend = "bottomright", 
                    ylegend = NULL, ...)
{
  if(inherits(ft, "DR"))
  {
    ft <- list(ft)
  }else if(!is.list(ft))
  {
    stop("argument ft must be a list of 'DR' objects")
  }else
  {
    if(any(sapply(ft, function(x) !inherits(x, "DR"))))        
      stop("argument ft must be a list of 'DR' objects")
  }
  
  nft <- length(ft)
  if (missing(datapch)) datapch <- 16
  if (missing(datacol)) datacol <- "black"
  if (missing(fitcol)) fitcol <- 2:(nft+1)
  if (missing(fitlty)) fitlty <- 1:nft
  fitcol <- rep(fitcol, length.out=nft)
  fitlty <- rep(fitlty, length.out=nft)
  
  if (missing(xlab))
    xlab <- "data"
  if (missing(ylab)) ylab <- "G(x)"
  if (missing(main)) main <- paste("Emp. and theo. exposure curve(s)")
  
  mydata <- ft[[1]]$data
  distname <- ft[[1]]$distname
  n <- length(mydata)
  s <- sort(mydata)
  largedata <- (n > 1e4)
  
  if(missing(xlim))
  {
    xmin <- min(mydata)
    xmax <- max(mydata)
    xlim <- c(xmin, xmax)
  }
  else
  {
    xmin <- xlim[1]
    xmax <- xlim[2]
  }
  
  verif.ftidata <- function(fti)
  {
    if (any(fti$data != mydata))
      stop("All compared fits must have been obtained with the same dataset")
    invisible()
  }
  lapply(ft, verif.ftidata)
  
  # computation of each fitted exposure curve
  sfin <- seq(xmin, xmax, length.out=101)
  comput.fti <- function(i)
  {
    fti <- ft[[i]]
    para <- c(as.list(fti$estimate), as.list(fti$fix.arg))
    distname <- fti$distname
    ecdistname <- paste("ec",distname,sep="")
    do.call(ecdistname, c(list(x=sfin), as.list(para)))
    
  }
  fittedec <- sapply(1:nft, comput.fti)
  
  #main plotting
  resec <- plot(eecf(x = mydata), main=main, xlab=xlab, ylab=ylab, xlim=xlim, 
                ylim=ylim, col=datacol, do.points=do.points)
  #plot fitted densities
  for(i in 1:nft)
    lines(sfin, fittedec[,i], lty=fitlty[i], col=fitcol[i], ...)
  
  if(addlegend)
  {
    if(missing(legendtext)) 
      legendtext <- paste("fit", 1:nft)
    legend(x=xlegend, y=ylegend, bty="n", legend=legendtext, 
           lty=fitlty, col=fitcol,...)
  }
}