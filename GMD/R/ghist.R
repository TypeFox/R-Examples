
gbreaks <-
  function(data,n)
{
  seq(min(data),max(data),length.out=n+1)
}


##' Generalized Histogram Computation with classes to contain a single histogram or multiple histograms
##'
##' \code{ghist} generates a single histogram.
##'
##' \code{gbreaks} generate bin boundaries for a histogram.
##'
##' \code{is.ghist} returns \code{TRUE} if \code{x} is an object of code{ghist} and
##' \code{FALSE} otherwise.
##' 
##' \code{as.ghist} is a generic function.  The method for numeric vectors will return a
##' \code{ghist} object.
##' 
##' \code{is.mhist} returns \code{TRUE} if \code{x} is an object of code{mhist} and
##' \code{FALSE} otherwise.
##' 
##' \code{as.mhist} is a generic function.  The method is for numeric list, matrices or data frames and will return a
##' \code{mhist} object.
##' 
##' \code{mhist2matrix} convert a \code{mhist} object into a numeric matrix, filling observations by row.
##'
##' 
##' @name ghist
##' @aliases ghist gbreaks is.ghist as.ghist mhist is.mhist as.mhist mhist2matrix
##' @title Generalized Histogram Computation
##' @usage
##' ghist(data,n,breaks=if (!.invalid(n)) NULL else "Sturges",
##' bins=NULL,digits=1)
##' 
##' gbreaks(data, n)
##' 
##' is.ghist(x)
##' 
##' as.ghist(x,bins)
##' 
##' is.mhist(x)
##' 
##' as.mhist(x,bins)
##'
##' mhist2matrix(h)
##' 
##' 
##' 
##' @param data a vector of values for which the histogram is desired.
##' @param n a single number giving the number of bins for the histogram.
##' @param breaks a vector giving the breakpoints between histogram bins, or
##' a character string naming an algorithm to compute the number of bins, or
##' a function to compute the number of bins (see \code{help("dist", package="graphics")}).
##' @param bins character vecter, the bin labels.
##' @param digits integer, the number of digits to round for breaks.
##' @param x an R object.
##' @param h an object of class \code{mhist}
##' @seealso \code{\link{plot.mhist}} \code{\link{mhist.summary}} \code{\link{plot.mhist.summary}}
##' @examples
##' ## load library
##' require("GMD")
##' 
##' ## create two normally-distributed samples
##' ## with unequal means and unequal variances
##' set.seed(2012)
##' v1 <- rnorm(1000,mean=-5, sd=10)
##' v2 <- rnorm(1000,mean=10, sd=5)
##' 
##' ## create common bins
##' n <- 20 # desired number of bins
##' breaks <- gbreaks(c(v1,v2),n) # bin boundaries
##' x <-
##'   list(ghist(v1,breaks=breaks,digits=0),
##'        ghist(v2,breaks=breaks,digits=0))
##' mhist.obj <- as.mhist(x)
ghist <-
  function(data,n,breaks=if (!.invalid(n)) NULL else "Sturges",bins=NULL,digits=1)
{
  if (!.invalid(n)) breaks <- gbreaks(data,n)
  
  res <- hist(data,breaks=breaks,plot=FALSE)
  ret <- res$counts
  breaks <- round(res$breaks,digits=digits)
  if (.invalid(bins)){
    ##attr(ret,"bins") <- paste("(",breaks[-length(breaks)],", ",breaks[-1],"]")
    attr(ret,"bins") <- breaks[-1]
  } else {
    attr(ret,"bins") <- bins
  }
  class(ret) <- c("ghist",class(ret))
  ret
}



is.ghist <-
  function(x)
{
  inherits(x,"ghist")
}


as.ghist <-
  function(x,bins)
{
  if(!is.numeric(x)){
    stop("`x' should be a numeric vector.")
  }
  if (.invalid(bins)){
    bins <- 1:length(x)
  }
  attr(x,"bins") <- bins
  class(x) <- c("ghist",class(x))
  x
}


is.mhist <-
  function(x)
{
  inherits(x,"mhist")
}


as.mhist <-
  function(x, bins=NULL)
{
  if(is.mhist(x)){
    return(x)
  } else if (inherits(x,"ghist")){
    x <- list(x)
  }
  
  if(!is.list(x) & !is.matrix(x) & !is.data.frame(x)){
    stop("`x' should be a list or matrix or data.frame; fail to turn it into a hist.")
  }

  if(is.list(x) & !is.data.frame(x)){
    x <- equalize.list(x)
    if(all(unlist(lapply(x,is.ghist)))){
      if (.invalid(bins)) bins <- attr(x[[1]],"bins")
      if (length(x)>1){
        if (!all(unlist(lapply(x,function(e) identical(bins,attr(e,"bins")))))){ # currently, bins should be in the same order.
          stop("To create a `mhist' (multi-histograms) object, each member histogram should has the same bins.")
        }
      }
    } else {
      if (.invalid(bins)) bins <- NULL
      n <- length(x[[1]])
    }
  } else {
    if (.invalid(bins)) bins <- colnames(x)
    if(is.matrix(x)) {
      x <- data.frame(x)
    }
    n <- ncol(x)
    tmp <- list()
    for (i in 1:nrow(x)){
      tmp[[i]] <- as.numeric(x[i,])
    }
    names(tmp) <- rownames(x)
    x <- tmp
  }

  if (is.null(bins)){
    bins <- 1:n
    warning("`bins' are not specified; use automatic numbering instead under assumption
that memeber histograms have the same order of the bins.")
  }
  
  if(!all(sapply(x, is.numeric))){
    stop("`x' should only contain numeric values; fail to turn it into a hist.")    
  }

  for (i in 1:length(x)) {
    attributes(x[[i]]) <- c()
  }
  class(x) <- c("mhist",class(x))
  attr(x,"memebers") <- length(x)
  if (!.invalid(names(x))){
    attr(x,"labels") <- names(x)
  } else {
    attr(x,"labels") <- as.character(1:length(x))
  }
  attr(x,"bins") <- bins
  x
}


mhist2matrix <-
  function(h)
{
  if (!is.mhist(h)){
    stop("`h' should be an object of class `mhist'.")
  }
  h <- equalize.list(h)
  ret <- do.call(rbind,h)
  rownames(ret) <- names(h)
  ret
}



##' Make member bins of a hist object equal size 
##' 
##' Make members of a list equal size 
##' @title Make members of a list equal size 
##' @param x a list of numeric vectors
equalize.list <-
  function(x)
{
  if (!is.list(x)){
    stop("`x' should be an object of class `list'.")
  }
  lv <- unlist(lapply(x,length))
  if(any(lv!=max(lv))){
    warning("Number of bins are not equal for all members. Shorter ones is padded with `NA' values.")
    x <- lapply(x,function(e) c(e,rep(NA,max(lv)-length(e))))
  }
  x
}



##' S3 method for class \code{mhist}
##'
##' Given a list, matrix or data.frame of histograms, plot multiple histograms side-by-side or as subplots.
##' @title S3 method for class `mhist'.
##' 
##' @export plot.mhist
##' 
##' @param x a numeric matrix or data frame, representing distributions by rows (bins by columns); or a list of numeric vectors as distributions. 
##' @param beside logical, whether plot histograms side-by-side.
##' @param labels a string vector of labels for the histograms in \code{x};
##' should have the same number as of the histograms.
##' @param colors the colors for the histograms; by default they are set to colors generated from palette Dark2.
##' Colors will be recycled if the size is smaller than the number of the histograms.
##' @param main an overall title for the plot. See \code{help("title", package="graphics")}.
##' @param sub a subtitle for the plot, describing the distance and/or alignment gap (the "shift").
##' @param ylab a title for the y axis. See \code{help("title", package="graphics")}.
##' @param xlab a title for the x axis. See \code{help("title", package="graphics")}.
##' @param xticks a string vector indicating the tickmark labels at x-axis. Defult: NULL.
##' @param xlabels character, labels at x-axis.
##' @param vlinePos numeric, posiitons for vertical lines.
##' @param x.las numeric in {0,1,2,3}; the style of axis labels.
##' See option \code{las} in\code{help("par", package="graphics")}.
##' @param xticks.type stinrg in {"pretty","original"}, whether plot the \code{xticks} in a pretty way or as is.
##' @param xlim range of x values, as in \code{help("plot", package="graphics")}.
##' @param ylim range of y values, as in \code{help("plot", package="graphics")}.
##' @param type type of plot, as in \code{help("plot", package="graphics")}.
##' @param font.type the name of a font type for drawing text. See \code{font} in \code{par}.
##' DEFAULT: \code{font.type = 1}, corresponding to plain text.
##' @param font.family the name of a font family for drawing text. See \code{family} in \code{par};
##' DEFAULT: \code{font.family = "sans"}, corresponding to san serif typeface.
##' @param cex.main a numerical value giving the amount by which \code{main}-title should be
##' magnified relative to the default.
##' @param cex.sub a numerical value giving the amount by which \code{sub}-title should be
##' magnified relative to the default.
##' @param cex.lab a numerical value giving the amount by which \code{xlab} and \code{ylab} should be
##' magnified relative to the default.
##' @param cex.tickmark a numerical value giving the amount by which tickmarks should be
##' magnified relative to the default.
##' @param cex.legend a numerical value giving the amount by which legends should be
##' magnified relative to the default.
##' @param tcl the length of tick marks as a fraction of the height of a line of text.
##' See option \code{tcl} in\code{help("par", package="graphics")}.
##' @param omi a vector of the form 'c(bottom, left, top, right)' giving the size of the outer margins in inches.
##' See option \code{omi} in\code{help("par", package="graphics")}.
##' @param mar a numerical vector of the form \code{c(bottom, left, top, right)} which gives the number of lines
##' of margin to be specified on the four sides of the plot.
##' See option \code{mar} in\code{help("par", package="graphics")}.
##' @param mgp the margin line (in 'mex' units) for the axis title, axis labels and axis line.
##' See option \code{mgp} in\code{help("par", package="graphics")}.
##' @param bin.unit numeric, indicating the width of a group of bar(s) in unit of x axis.
##' @param legend.lab legend labels, a string vector of the same length of distributions in \code{x},
##' using \code{labels} by default. No legend is displayed when it is \code{NA}. 
##' @param legend.pos string, a keyword to be used to position the legend.
##' See \code{help("legend", package="graphics")}.
##' @param ... arguments to be passed to method \code{plot.mhist},
##' such as graphical parameters (see \code{par}).
##' @references See \code{help(GMD)}
##' @seealso \code{\link{mhist}} \code{\link{mhist.summary}} \code{\link{plot.mhist.summary}} \code{\link{plot.gmdp}} \code{\link{plot.gmdm}}
##' @keywords methods hplot
##' @examples
##' ## load library
##' require("GMD")
##' 
##' ## create two normally-distributed samples
##' ## with unequal means and unequal variances
##' set.seed(2012)
##' v1 <- rnorm(1000,mean=-5, sd=10)
##' v2 <- rnorm(1000,mean=10, sd=5)
##' 
##' ## create common bins
##' n <- 20 # desired number of bins
##' breaks <- gbreaks(c(v1,v2),n) # bin boundaries
##' x <-
##'   list(ghist(v1,breaks=breaks,digits=0),
##'        ghist(v2,breaks=breaks,digits=0))
##' mhist.obj <- as.mhist(x)
##' 
##' ## plot histograms side-by-side
##' plot(mhist.obj,mar=c(1.5,1,1,0),
##' main="Histograms of simulated normal distributions")
##' 
##' ## plot histograms as subplots,
##' ## with corresponding bins aligned
##' plot(mhist.obj,beside=FALSE,mar=c(1.5,1,1,0),
##'      main="Histograms of simulated normal distributions")
plot.mhist <-
  function(x,
           ##
           beside=TRUE,
           ## labels and color
           labels=NULL,
           colors=NULL,
           ## title and axis-labels
           main=NULL,
           sub=NULL,
           ylab=NULL,
           xlab=NULL,
           xticks=NULL,
           xlabels=NULL,
           vlinePos=NULL,
           x.las=1,
           xticks.type=c("pretty","original"),
           xlim=NULL,
           ylim=NULL,
           type=NULL, ##XB
           ## font and size
           font.type=1,
           font.family=c("sans","serif","mono"),           
           cex.main=1.75,         
           cex.sub=cex.main*0.9,
           cex.lab=1.25,
           cex.tickmark=0.75,
           cex.legend=1.5,
           ## par
           tcl=-0.25,
           omi=c(0.5,0.5,1.0,0.25),
           mar=c(4,1,0,1),
           mgp=c(0,0.5,0),
           ## graphical param of the bin
           bin.unit=0.8,
           ## legend
           legend.lab=labels, 
           legend.pos=c("topright","top","topleft"),
           ...
           )
{
  message(xlabels)
  
  if (!is.mhist(x)){
    stop("`x' should be an object of class `mhist'.")
  }
  
  ##
  legend.pos <- match.arg(legend.pos)
  font.family <- match.arg(font.family)
  xticks.type <- match.arg(xticks.type)
  half.bin.unit <- bin.unit/2
  

  ##
  ret <- list()

  x.ori <- x
  bins <- attr(x,"bins")
  x <- do.call(cbind, x)

  ##
  M <- nrow(x) # number of bins
  N <- ncol(x) # number of histograms

  ## par
  ## labels
  if (.invalid(labels)){
    labels <- colnames(x)
  }
  if (is.null(labels)){
    labels <- 1:N
  }

  
  ## colors
  if (.invalid(colors)){
    colors <- .brewer.pal.Dark2(8)
  }
  
  if(length(colors) < N){
    warning("The length of `colors' is shorter than that of `x' and the `colors' are recycled. Specifying a longer vector of colors is recommended!")
    colors <- rep(colors,ceiling(N/length(colors)))
  }

  ## main and sub ##
  if (.invalid(main)){
    if(length(x.ori) == 1)
      main <- "Histogram"
    else
      main <- "Multiple Histograms"
  }
  
  if (.invalid(sub)){
    sub <- ""
  }

  
  ## xlim ##
  if (.invalid(xlim))
    xlim <- c(1, M)
  xlim.extra <- c(xlim[1]-(1-half.bin.unit), xlim[2]+(1-half.bin.unit))
  
  ## ylim ##
  if (.invalid(ylim))
    ylim <- .scale.range(range(unlist(x), na.rm=TRUE),1.05)
  y.neg <- any(x<0, na.rm=TRUE)

  ##
  if (.invalid(xlab)){
    xlab <- "Bin"
  }

  ##
  if (.invalid(ylab)){
    ylab <- "Count"
  }

  ## at ##
  binsize <- length(bins)
  x.at <- 1:binsize
  
  if (is.null(xticks)){
    if (is.numeric(bins)){
      if (xticks.type=="pretty"){
        x.which <- pretty(c(1,binsize))
        x.labels <- bins[x.which]
        x.at <- x.at[x.which]
      }
    } else {
      x.at <- x.at
      x.labels <- bins
    }
  } else if (identical(xticks,NA)) {
    x.at <- x.labels <- NULL
  } else if (is.numeric(xticks)) {
    x.at <- xticks
    if (.invalid(xlabels)){
      x.labels <- NULL
    } else{
      x.labels <- xlabels
    }
  } 

  ##print(x.at)
  y.at <- pretty(ylim)

  ##par
  ##old.par <- par(no.readonly=TRUE)
  ##on.exit(par(old.par))

  ## par(font.main=font.type, font.lab=font.type, font.axis=font.type)
  ## par(family=font.family)
  ## par(omi=omi,mar=mar,mgp=mgp,tcl=tcl)
  
  if (beside){
    plot(xlim.extra,
         ylim,
         type="n",
         main="",
         xlab="",
         ylab="",
         xlim=xlim.extra,
         ylim=ylim,
         xaxt="n",
         yaxt="n",
         cex.lab=cex.lab,
         ...
         )

    if (!is.null(vlinePos)){
      abline(v=vlinePos,lty="dotdash")
    }
    
    mtext(line=2,side=3,text=main,outer=FALSE,cex=cex.main,adj=0) # main title
    mtext(line=0.55,side=3,text=sub,outer=FALSE,cex=cex.sub,adj=0) # subtitle

    ##print(x.at)
    if (!is.null(x.at))
      axis(side=1, at=x.at, labels=x.labels, cex.axis=cex.tickmark,las=x.las)
    
    if (!is.null(y.at))
      axis(2, at=y.at, cex.axis=cex.tickmark, las=1)

    bottoms <- 0
    bottoms <- matrix(bottoms, nrow=M, ncol=N)
    barhalf.bin.unit <- 2*half.bin.unit/N
    border <- NA

    for (i in 1:N){
      if(!length(type)){
        rect(xleft=(1:M)-half.bin.unit+(i-1)*barhalf.bin.unit,
             ybottom=bottoms[,i],
             xright=(1:M)-half.bin.unit+i*barhalf.bin.unit,
             ytop=x[,i],
             col=colors[i],
             border=border)
      } else if (type %in% "polygon"){
        .x <- (1:M)-half.bin.unit+(2*i-1)/2*barhalf.bin.unit
        .y <- x[,i]
        y.polygon <- c(.y, rep(0,length(.y)))
        x.polygon <- c(.x, rev(.x))
        graphics::polygon(x.polygon, y.polygon, col=colors[i], border=border)        
      } else{
        lines(.x,x[,i],col=colors[i],type=type)
      }
    }
    
    if (!all(.invalid(legend.lab))) {
      ##print("Plotting legends ..")
      xjust <- yjust <- 0.5
      legend(legend.pos,
             legend=legend.lab,
             text.col="black",
             cex=cex.legend,
             box.col=FALSE,
             fill=colors,
             border=border,
             xjust=xjust,
             yjust=yjust)
    }
    ## xlab and ylab
    mtext(line=3,side=1,text=xlab,outer=FALSE,cex=cex.lab,adj=0.5)
    mtext(line=3,side=2,text=ylab,outer=FALSE,cex=cex.lab,adj=0.5) 

  } else {
    par(mfrow=c(N,1))
    par(omi=omi,mar=mar,mgp=mgp,tcl=tcl)
    for (j in 1:N){
      i.x <- x.ori[[j]]
      plot(xlim.extra,
           ylim,
           type="n",
           main="",
           xlab="",
           ylab="",
           xlim=xlim.extra,
           ylim=ylim,
           xaxt="n",
           yaxt="n",
           cex.lab=cex.lab,
           ...
           )
      
      if (!is.null(vlinePos)){
        abline(v=vlinePos,lty="dotdash")
      }
      ##print(x.at)
      if (!is.null(x.at))
        axis(side=1, at=x.at, labels=x.labels, cex.axis=cex.tickmark,las=x.las)
      
      if (!is.null(y.at))
        axis(2, at=y.at, cex.axis=cex.tickmark, las=1)
      
      bottoms <- 0
      bottoms <- matrix(bottoms, nrow=M, ncol=1)
      barhalf.bin.unit <- 2*half.bin.unit/1
      border <- NA
      i <- 1
      rect(xleft=(1:M)-half.bin.unit+(i-1)*barhalf.bin.unit,
           ybottom=bottoms[,i],
           xright=(1:M)-half.bin.unit+i*barhalf.bin.unit,
           ytop=i.x,
           col=colors[j],
           border=border) 

      if (!all(.invalid(legend.lab))) {
        ##print("Plotting legends ..")
        xjust <- yjust <- 0.5
        legend(legend.pos,
               legend=legend.lab[j],
               text.col="black",
               cex=cex.legend,
               box.col=FALSE,
               fill=colors[j],
               border=border,
               xjust=xjust,
               yjust=yjust)
      }
    }
    
    mtext(line=2,side=3,text=main,outer=TRUE,cex=cex.main,adj=0) # main title
    mtext(line=0.55,side=3,text=sub,outer=TRUE,cex=cex.sub,adj=0) # subtitle
    ## xlab and ylab
    mtext(line=1,side=1,text=xlab,outer=TRUE,cex=cex.lab,adj=0.5)
    mtext(line=1,side=2,text=ylab,outer=TRUE,cex=cex.lab,adj=0.5) 
  }

  
  ret$members <- N
  ret$bins <- M
  ret$xlim <- xlim
  ret$xlim.extra <- xlim.extra
  ret$ylim <- ylim
  ret$x.at <- x.at
  ret$y.at <- y.at

  invisible(ret)
}



##' Bin-wise summary of a \code{mhist} object of histograms
##'
##' Bin-wise summary of a \code{mhist} object of histograms
##' @name mhist.summary 
##' @aliases mhist.summary plot.mhist.summary
##' 
##' @title Bin-wise summary of histograms
##' 
##' @usage
##' mhist.summary(h, ...)
##' 
##' \method{plot}{mhist.summary}(x,bins,plot.ci=TRUE,col=NULL,
##' ci.color="orchid1",tcl=-0.25,omi=c(0.5,0.5,1.0,0.25),mar=c(3,3,3,1),
##' mgp=c(2,0.5,0),if.plot.new=TRUE,...)
##' 
##' @param h a \code{"mhist"} object as produced by \code{as.mhist}
##' @param x a \code{mhist.summary} object as produced by \code{mhist.summary}
##' @param bins character vecter, the bin labels; if non-specific, bins are numbered/labeled starting with one.
##' @param plot.ci logical, indicating whether plot error bars that represent the 0.50 confidence interval (CI) 
##' @param col color of the histogram
##' @param ci.color color of the error bars
##' @param tcl the length of tick marks as a fraction of the height of a line of text.
##' See option \code{tcl} in\code{help("par", package="graphics")}.
##' @param omi a vector of the form 'c(bottom, left, top, right)' giving the size of the outer margins in inches.
##' See option \code{omi} in\code{help("par", package="graphics")}.
##' @param mar a numerical vector of the form \code{c(bottom, left, top, right)} which gives the number of lines
##' of margin to be specified on the four sides of the plot.
##' See option \code{mar} in\code{help("par", package="graphics")}.
##' @param mgp the margin line (in 'mex' units) for the axis title, axis labels and axis line.
##' @param if.plot.new logical, whether starting a new device or not.
##' @param ... arguments to be passed to method \code{plot.mhist.summary}.\cr
##' See \code{help("barplot2", package="gplots")}.
##' @return
##' \code{mhist.summary} returns a \code{mhist.summary} object
##' @seealso \code{\link{mhist}} \code{\link{plot.mhist}} \code{\link{plot.gmdp}} \code{\link{plot.gmdm}}
mhist.summary <-
  function(h, ...)
{
  if (!is.mhist(h)){
    stop("`h' should be an object of class `mhist'.")
  }
  bins <- attr(h,"bins")
  h <- mhist2matrix(h)

  ret <- matrix(nrow=7,ncol=ncol(h)) ## features-by-bins
  rownames(ret) <- c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","NA's")
  colnames(ret) <- c()
  tmp.apply <- apply(h, 2, function(e,...) summary(e,...))
  
  if (is.list(tmp.apply)){
    for (i in 1:length(tmp.apply)){
      e <- tmp.apply[[i]]
      ret[names(e),i] <- e
    }
  } else {
    ret <- tmp.apply
  }
  class(ret) <- c("mhist.summary",class(ret))
  attr(ret,"members") <- nrow(h)
  attr(ret,"bins") <- bins
  attr(ret,"nbin") <- ncol(h)
  ret
}



plot.mhist.summary <-
  function(x,
           bins,
           plot.ci=TRUE,
           col=NULL,
           ci.color="orchid1",
           ## graphical params
           tcl=-0.25,
           omi=c(0.5,0.5,1.0,0.25),
           mar=c(3,3,3,1),
           mgp=c(2,0.5,0),
           if.plot.new=TRUE,
           ...)
{
  if (!inherits(x,"mhist.summary")){
    stop("`x' should be an object of class `mhist.summary'.")
  }
  if (.invalid(bins)) bins <- attr(x,"bins")
  if (if.plot.new) {
    par(omi=omi,mar=mar,mgp=mgp,tcl=tcl)
  }
  ##print(bins)
  barplot2(x["Mean",],
           names.arg=bins,
           ci.l=x["1st Qu.",],
           ci.u=x["3rd Qu.",],
           plot.ci=plot.ci,
           col=col,
           ci.color=ci.color,
           ...)
}

