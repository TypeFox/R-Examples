# $Id: ooplot.R 1557 2012-06-08 17:56:37Z warnes $

ooplot <- function(data, ...) UseMethod("ooplot")

ooplot.default <- function(data, width=1, space=NULL, names.arg=NULL,
                           legend.text=NULL, horiz=FALSE,
                           density=NULL, angle=45, kmg="fpnumkMGTP",
                           kmglim=TRUE,
                           type=c("xyplot", "linear", "barplot", "stackbar"),
                           col=heat.colors(NC), prcol=NULL,
                           border=par("fg"), main=NULL, sub=NULL,
                           xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL,
                           xpd=TRUE, log="", axes=TRUE,
                           axisnames=TRUE, prval=TRUE, lm=FALSE,
                           cex.axis=par("cex.axis"),
                           cex.names=par("cex.axis"),
                           cex.values=par("cex"),inside=TRUE,
                           plot=TRUE, axis.lty=0, plot.ci=FALSE,
                           ci.l=NULL, ci.u=NULL, ci.color="black",
                           ci.lty="solid", ci.lwd=1, plot.grid=FALSE,
                           grid.inc=NULL, grid.lty="dotted",
                           grid.lwd=1, grid.col="black", add=FALSE,
                           by.row=FALSE, ...)
{

  ##
  ##  oopplot function block
  ##
  ## this is the location of the helper functions for this method
  ##


  optlim <- function(lim, log=FALSE) {
    ## define xlim and ylim, adjusting for log-scale if needed
    factor <- 1.05
    if (log) {
      min <- 10^floor(log10(lim[1]/factor))
      max <- 10^ceiling(log10(factor*lim[2]))
    } else {
      range <- (factor*lim[2]-lim[1]/factor)
      if (range>0) {
        ## we know the range, now find the optimal start and endpoints
        scale <- 10^floor(log10(range))
        min <- scale*(floor((lim[1]/factor)/scale))
        max <- scale*(ceiling(factor*lim[2]/scale))
      } else {
        min=0
        max=1
      }
    }
    if (type=="barplot" || type=="stackbar") max=max+1

    return(c(min, max))
  }

  linearfitplot <- function(x,y,xlim,col) {
    ## calculate a linear fit through the datapoints and plot
    local <- data.frame(x=x,y=y)
    local.lm <- lm(y ~ x, data=local)
    summary <- summary(local.lm)
    xmin=min(x)
    xmax=max(x)
    if (xlim[1]<xmin) {
      x1=mean(xlim[1],xmin)
    } else {
      x1=max(min(x),xlim[1])
    }
    if (xlim[2]>xmax) {
      x2=mean(xlim[2],xmax)
    } else {
      x2=min(max(x),xlim[2])
    }
    y1=summary$coefficients[2,1]*x1+summary$coefficients[1,1]
    y2=summary$coefficients[2,1]*x2+summary$coefficients[1,1]
    p1=c(x1,x2)
    p2=c(y1,y2)
    lines(p1,p2,col=col,lty=2,lwd=2)
  }


  castNA <- function(matrix) {
    newmatrix <- matrix
    for (j in 1:ncol(matrix)) {
      for (i in 1:nrow(matrix)) {
        newmatrix[i,j] <- ifelse(is.na(matrix[i,j]),0,matrix[i,j])
      }
    }
    return(newmatrix)
  }
  ##
  ## End of function block
  ##
  ## In R, most people think about data in columns rather than rows
  if(by.row)
    data <- as.matrix(data)
  else
    data <- t(as.matrix(data))

  ## make sure we only accept the supported plot options
  type <- match.arg(type)
  ## check data validity
  if( nrow(data) < 2 ) stop("At least 2 columns are required.")
  ## check defaults
  if (!missing(inside)) .NotYetUsed("inside", error=FALSE)# -> help(.)
  ## set the beside parameter
  if (type=="stackbar") {
    beside <- FALSE
    data <- castNA(data)
  }
  else
    beside <- TRUE

  ## split the data into x and y values
  height <- data[-1,,drop=FALSE]


  heightscale <- ""
  heightsymbol <- ""
  ##
  if ((kmg!="") && (kmglim==TRUE)){
    ##
    ## auto scale the parameters
    ##
    ## now scale the data, valid factors are
    ## P : peta=1E15
    ## T : tera=1E12
    ## G : giga=1E09
    ## M : mega=1E06
    ## k : kilo=1E03
    ## m : milli=1E-03
    ## u : micro=1E-06
    ## n : nano=1E-09
    ## p : pico=1E-12
    ## f : femto=1E-15
    maxheight <- max(abs(height), na.rm=TRUE)
    heightfactor <- 1
    if (maxheight>10E15 && any(grep("P", kmg))) {
      heightfactor <- 1E15
      heightscale <- "Peta"
      heightsymbol <- "P"
    }
    else if (maxheight>1E12 && any(grep("T", kmg))) {
      heightfactor <- 1E12
      heightscale <- "Tera"
      heightsymbol <- "T"
    }
    else if (maxheight>1E09 && any(grep("G", kmg))) {
      heightfactor <- 1E09
      heightscale <- "Giga"
      heightsymbol <- "G"
    }
    else if (maxheight>1E06 && any(grep("M", kmg))) {
      heightfactor <- 1E06
      heightscale <- "Mega"
      heightsymbol <- "M"
    }
    else if (maxheight>1E03 && any(grep("k", kmg))) {
      heightfactor <- 1E03
      heightscale <- "Kilo"
      heightsymbol <- "k"
    }
    else if (maxheight<1E-15 && any(grep("f", kmg))) {
      heightfactor <- 1E-15
      heightscale <- "Femto"
      heightsymbol <- "f"
    }
    else if (maxheight<1E-12 && any(grep("p", kmg))) {
      heightfactor <- 1E-12
      heightscale <- "Pico"
      heightsymbol <- "p"
    }
    else if (maxheight<1E-09 && any(grep("n", kmg))) {
      heightfactor <- 1E-09
      heightscale <- "Nano"
      heightsymbol <- "n"
    }
    else if (maxheight<1E-06 && any(grep("u", kmg))) {
      heightfactor <- 1E-06
      heightscale <- "Micro"
      heightsymbol <- "u"
    }
    else if (maxheight<1E-03 && any(grep("m", kmg))) {
      heightfactor <- 1E-03
      heightscale <- "Milli"
      heightsymbol <- "m"
    }

    height <- height/heightfactor

  }

  ## fill the xaxis data set and matching rownames
  xaxis <- data[1, ]
  rownames <- rownames(data)

  ##
  if (missing(space))
    space <- if (is.matrix(height) && beside) c(0, 1) else 0.2
  space <- space * mean(width)
  ##
  if (plot && axisnames && missing(names.arg)) {
    if (type=="xyplot")
      names.arg <- colnames(height)
    else if (type=="linear")
      names.arg <- xaxis
    else if (type=="barplot")
      names.arg <- xaxis
    else if (type=="stackbar")
      names.arg <- xaxis
    else names.arg <- xaxis
  }


  ## set the legend text if it is null
  if(is.logical(legend.text))
    legend.text <- rownames[-c(1)]

  ##  if(legend.text && is.matrix(height)) rownames(height)
  ##  else colnames(height)

  if (is.vector(height) || is.array(height))
    {
      height <- rbind(height)
    }
  else if (!is.matrix(height))
    stop("`height' must be a vector or a matrix")

  ## Check for log scales
  logx <- FALSE
  logy <- FALSE

  if (log !="")
    {
      if (any(grep("x", log)))
        logx <- TRUE
      if (any(grep("y", log)))
        logy <- TRUE
    }

  ## Cannot "hatch" with rect() when log scales used
  if ((logx || logy) && !is.null(density))
    stop("Cannot use shading lines in bars when log scale is used")

  ## Set the size of Rows and Columns
  NR <- nrow(height)
  NC <- ncol(height)
  ## w.r, w.l, w.m are the x-axis coordinates

  if (type=="barplot") {
    if (NR<1) NR <- 1
    if (NC<1) NC <- 1
    if (length(space)==2) {
      space <- rep.int(c(space[2], rep.int(space[1], NR - 1)), NC)
    }
    width <- rep(width, length=NR * NC)
  } else if (type=="stackbar") {
    width <- rep(width, length=NC)
  }  else {
    width <- rep(width, length=NC)
    delta <- width / 2
  }

  ## set the proper x-axis scale
  ## linear is a switch between equidistant and scaled
  if (type=="xyplot") {
    delta <- 0
    w.m <- xaxis
    w.r <- w.m + delta
    w.l <- w.m - delta
  }
  else if (type=="linear") {
    w.r <- cumsum(width)
    w.m <- w.r - delta
    w.l <- w.m - delta
    xaxis <- w.m
  }
  else if(type=="barplot" || type=="stackbar") {
    delta <- width / 2
    w.r <- cumsum(space + width)
    w.m <- w.r - delta
    w.l <- w.m - delta
    ##if graphic will be stacked bars, do not plot ci
    if (beside && (NR > 1) && plot.ci)
      plot.ci=FALSE
  }
  else stop("Unkown plot type")


  ## error check ci arguments
  if (plot && plot.ci)
    {
      if ((missing(ci.l)) || (missing(ci.u)))
        stop("confidence interval values are missing")

      if (is.vector(ci.l))
        ci.l <- cbind(ci.l)
      else if (is.array(ci.l) && (length(dim(ci.l))==1))
        ci.l <- rbind(ci.l)
      else if (!is.matrix(ci.l))
        stop("`ci.l' must be a vector or a matrix")

      if (is.vector(ci.u))
        ci.u <- cbind(ci.u)
      else if (is.array(ci.u) && (length(dim(ci.u))==1))
        ci.u <- rbind(ci.u)
      else if (!is.matrix(ci.u))
        stop("`ci.u' must be a vector or a matrix")

      if ( any(dim(height) !=dim(ci.u) ) )
        stop("'height' and 'ci.u' must have the same dimensions.")
      else if ( any( dim(height) !=dim(ci.l) ) )
        stop("'height' and 'ci.l' must have the same dimensions.")
    }

  if (beside)
    w.m <- matrix(w.m, ncol=NC)

  ## check height/ci.l if using log scale to prevent log(<=0) error
  ## adjust appropriate ranges and bar base values
  if ((logx && horiz) || (logy && !horiz))
    {
      if (min(height, na.rm=TRUE) <=0)
        stop("log scale error: at least one 'height' value <=0")

      if (plot.ci && (min(ci.l) <=0))
        stop("log scale error: at least one lower c.i. value <=0")

      if (logx && !is.null(xlim) && (xlim[1] <=0))
        stop("log scale error: 'xlim[1]' <=0")

      if (logy && !is.null(ylim) && (ylim[1] <=0))
        stop("'log scale error: 'ylim[1]' <=0")

      ## arbitrary adjustment to display some of bar for min(height)
      ## or min(ci.l)
      if (plot.ci)
        rectbase <- rangeadj <- (0.9 * min(c(height, ci.l)))
      else
        rectbase <- rangeadj <- (0.9 * min(height))

      ## if axis limit is set to < above, adjust bar base value
      ## to draw a full bar
      if (logy && !is.null(ylim))
        rectbase <- ylim[1]
      else if (logx && !is.null(xlim))
        rectbase <- xlim[1]

      ## if stacked bar, set up base/cumsum levels, adjusting for log scale
      if (type=="stackbar") {
        heightdata <- height ## remember the original values
        height <- rbind(rectbase, apply(height, 2, cumsum))
      }

      ## if plot.ci, be sure that appropriate axis limits are set to
      ## include range(ci)
      lim <-
        if (plot.ci)
          c(height, ci.l, ci.u)
        else
          height
    }
  else
    {
      ## Use original bar base value
      rectbase <- 0

      ## if stacked bar, set up base/cumsum levels
      if (type=="stackbar") {
        heightdata <- height ## remember the original values

        height <- rbind(rectbase, apply(height, 2, cumsum))
      }

      ## if plot.ci, be sure that appropriate axis limits are set to
      ## include range(ci)
      lim <-
        if (plot.ci)
          c(height, ci.l, ci.u)
        else
          height

      ## use original range adjustment factor
      rangeadj <- (-0.01 * lim)
    }


  ## calculate the ranges ourselves
  if (missing(xlim)) xlim <- optlim(range(w.l, w.r, na.rm=TRUE), logx)
  if (missing(ylim)) ylim <- optlim(range(height, na.rm=TRUE), logy)
  ##
  if(plot) ##-------- Plotting :
    {
      if (type=="barplot" || type=="stackbar") {
        if (horiz)

          rectbase <- xlim[1] ## make sure the xlimit and the rectbase
                              ## are in sync
        else
          rectbase <- ylim[1] ## make sure the ylimit and rectbase are in sync
        opar <-
          if (horiz)
            par(xaxs="i", xpd=xpd)
          else
            par(yaxs="i", xpd=xpd)
      } else {
        opar <-
          if (horiz)
            par(xaxs="r", yaxs="r", xpd=xpd)
          else
            par(xaxs="r", yaxs="r", xpd=xpd)
      }
      on.exit(par(opar))

      ## If add=FALSE open new plot window
      ## else allow for adding new plot to existing window
      if (!add)
        {
          plot.new()
          plot.window(xlim, ylim, log=log, ...)
        }

      ## Set plot region coordinates
      usr <- par("usr")

      ## adjust par("usr") values if log scale(s) used
      if (logx)
        {
          usr[1] <- 10 ^ usr[1]
          usr[2] <- 10 ^ usr[2]
        }

      if (logy)
        {
          usr[3] <- 10 ^ usr[3]
          usr[4] <- 10 ^ usr[4]
        }

      ## if prcol specified, set plot region color
      if (!missing(prcol))
        rect(usr[1], usr[3], usr[2], usr[4], col=prcol)

      ## if plot.grid, draw major y-axis lines if vertical or x axis
      ## if horizontal R V1.6.0 provided axTicks() as an R equivalent
      ## of the C code for CreateAtVector.  Use this to determine
      ## default axis tick marks when log scale used to be consistent
      ## when no grid is plotted.  Otherwise if grid.inc is specified,
      ## use pretty()

      if (plot.grid)
      {
        par(xpd=FALSE)

        if (is.null(grid.inc))
        {
          if (horiz)
          {
            grid <- axTicks(1)
            abline(v=grid, lty=grid.lty, lwd=grid.lwd, col=grid.col)
          }
          else
          {
            grid <- axTicks(2)
            abline(h=grid, lty=grid.lty, lwd=grid.lwd, col=grid.col)
          }
        }
        else
        {
          if (horiz)
          {
            grid <- pretty(xlim, n=grid.inc)
            abline(v=grid, lty=grid.lty, lwd=grid.lwd, col=grid.col)
          }
          else
          {
            grid <- pretty(ylim, n=grid.inc)
            abline(h=grid, lty=grid.lty, lwd=grid.lwd, col=grid.col)
          }
        }

         par(xpd=xpd)
      }

      ##
      ## end of the general setup, now get ready to plot the elements
      ##
      ## cycle through all the sets and plot the lines
      if (type=="xyplot" || type=="linear") {
        pch <- c()
        for (i in 1:NR) {
          list <- is.finite(height[i, ])
          xrange <- xaxis[list]
          yrange <- height[i, list]
          lines(xrange, yrange, col=col[i])
          if ((type=="xyplot") && (lm==TRUE)) linearfitplot(xrange,yrange,xlim,col[i])
          symbol=21+(i %% 5)
          points(xrange, yrange, pch=symbol, bg=col[i], col=col[i])
          pch <- c(pch, symbol)
          ##
          if (prval) {
            values=paste(format(yrange, digits=4), heightsymbol, sep="")
            chh <- par()$cxy[2]
            chw <- par()$cxy[1]
            if (logy) {
              factor=(yrange+chh)/yrange
              text(xrange, 1.1*yrange, labels=values, adj=c(0, 1),
	           srt=0, cex=cex.values, bg="white", col="navy")
              } else {
              text(xrange, yrange+chh, labels=values, adj=c(0, 1),
	           srt=0, cex=cex.values, bg="white", col="navy")
            }
          }
        }
      }

      if (type=="barplot" || type=="stackbar") {
        xyrect <- function(x1, y1, x2, y2, horizontal=TRUE, ...)
          {
            if(horizontal)
              rect(x1, y1, x2, y2, ...)
            else
              rect(y1, x1, y2, x2, ...)
          }

        chh <- par()$cxy[2]
        chw <- par()$cxy[1]
        if (type=="barplot") {
          xyrect(rectbase, w.l, c(height), w.r, horizontal=horiz,
                 angle=angle, density=density, col=col, border=border)
          if (prval) {
            values=paste(format(c(height), digits=4), heightsymbol, sep="")
                if (horiz) {
                   text(c(height)+chw, w.l, labels=values, adj=c(0, 1),
                        srt=0, cex=cex.values, bg="white", col="navy")
                } else {
                  text(w.l, c(height)+chh, labels=values, adj=c(0, 1),
                       srt=0, cex=cex.values, bg="white", col="navy")
         }
          }
        } else if (type=="stackbar") {
          for (i in 1:NC) {
            xyrect(height[1:NR, i], w.l[i], height[-1, i], w.r[i],
                   horizontal=horiz, angle=angle, density=density,
                   col=col, border=border)
            if (prval) {
              for (j in 1:NR) {
                values=paste(format(c(heightdata[j, i]), digits=4),
                  heightsymbol, sep="")
                if (horiz) {
                   text(c(height[j, i])+chw, w.m[i], labels=values,
                        adj=c(0, 1), srt=0, cex=cex.values, bg="white",
                        col="navy")
                } else {
                text(w.m[i], c(height[j, i])+chh, labels=values, adj=c(0, 1),
                     srt=0, cex=cex.values, bg="white", col="navy")
              }
              }
            }
          }
        }
      }

      if (plot.ci)
        {
          ## CI plot width=barwidth / 2
          ci.width=width / 4

          if (horiz)
            {
              segments(ci.l, w.m, ci.u, w.m, col=ci.color, lty=ci.lty,
                       lwd=ci.lwd)
              segments(ci.l, w.m - ci.width, ci.l, w.m + ci.width,
                       col=ci.color, lty=ci.lty, lwd=ci.lwd)
              segments(ci.u, w.m - ci.width, ci.u, w.m + ci.width,
                       col=ci.color, lty=ci.lty, lwd=ci.lwd)
            }
          else
            {
              segments(w.m, ci.l, w.m, ci.u, col=ci.color, lty=ci.lty,
                       lwd=ci.lwd)
              segments(w.m - ci.width, ci.l, w.m + ci.width, ci.l,
                       col=ci.color, lty=ci.lty, lwd=ci.lwd)
              segments(w.m - ci.width, ci.u, w.m + ci.width, ci.u,
                       col=ci.color, lty=ci.lty, lwd=ci.lwd)
            }
        }


      if (axisnames && !is.null(names.arg)) # specified or from {col}names
        {
          at.l <-
            if (length(names.arg) !=length(w.m))
              {
                if (length(names.arg)==NC) # i.e. beside (!)
                  colMeans(w.m)
                else
                  if ((type=="barplot") && (NR==1)) {
                    median(w.m)
                  } else if ((type=="linear") && (NR==1)) {
                    median(w.m)
                  } else {
                    stop("incorrect number of names now")
                  }
              }
            else w.m

          ##axis(if(horiz) 2 else 1, at=at.l, labels=names.arg,
          ##     lty=axis.lty, cex.axis=cex.names, ...)
        }


      if(!is.null(legend.text))
        {
          legend.col <- rep(col, length=length(legend.text))

          if((horiz & beside) || (!horiz & !beside))
            {
              legend.text <- rev(legend.text)
              legend.col <- rev(legend.col)
              density <- rev(density)
              angle <- rev(angle)
            }

          ## adjust legend x and y values if log scaling in use
          if (logx)
            legx <- usr[2] - ((usr[2] - usr[1]) / 10)
          else
            legx <- usr[2] - xinch(0.1)

          if (logy)
            legy <- usr[4] - ((usr[4] - usr[3]) / 10)
          else
            legy <- usr[4] - yinch(0.1)

          if (type=="barplot" || type=="stackbar") {
            legend(legx, legy,
                   legend=legend.text, angle=angle, density=density,
                   fill=legend.col, xjust=1, yjust=1)
          } else {
            legend(legx, legy,
                   legend=legend.text,
                   ##angle=angle,
                   ##density=density,
                   lwd=1,
                   col=legend.col,
                   pch=pch,
                   pt.bg=legend.col,
                   xjust=1,
                   yjust=1)
          }
        }

      title(main=main, sub=sub, xlab=xlab, ylab=paste(heightscale, ylab), ...)

      ## if axis is to be plotted, adjust for grid "at" values
      if (axes)
        {
          par(lab=c(10, 10, 7))
          if (type=="barplot" || type=="stackbar") {
            axis(if(horiz) 2 else 1, at=at.l, labels=names.arg, lty=axis.lty,
                 cex.axis=cex.names, ...)
            if(plot.grid)
              axis(if(horiz) 1 else 2, at=grid, cex.axis=cex.axis, ...)
            else
              axis(if(horiz) 1 else 2, cex.axis=cex.axis, ...)
          }
          else if (type=="linear") {
            if(plot.grid) {
              axis(if(horiz) 1 else 2, at=grid, cex.axis=cex.axis, ...)
              axis(if(horiz) 2 else 1, at=at.l, labels=names.arg,
                   lty=axis.lty, cex.axis=cex.names, ...)
            } else {
              axis(if(horiz) 1 else 2, pos=0, cex.axis=cex.axis, ...)
              axis(if(horiz) 2 else 1, pos=ylim[1], at=at.l, labels=names.arg,
                   cex.axis=cex.names, ...)
            }
          }
          else if (type=="xyplot") {
            if(plot.grid) {
              axis(if(horiz) 1 else 2, at=grid, cex.axis=cex.axis, ...)
              axis(if(horiz) 2 else 1, cex.axis=cex.axis, )
            } else {
              if (horiz) {
                axis(1, pos=ylim[1], cex.axis=cex.axis, ...)
                axis(2, pos=xlim[1], cex.axis=cex.axis, ...)
              } else {
                axis(2, pos=xlim[1], cex.axis=cex.axis, ...)
                axis(1, pos=ylim[1], cex.axis=cex.axis, ...)
              }
            }
          }
        }

      invisible(w.m)

    }

    else w.m
}
