# $Id: balloonplot.R 1302 2008-08-14 20:06:22Z warnes $

balloonplot <- function(x,...)
  UseMethod("balloonplot",x)

balloonplot.table <- function(x, xlab, ylab, zlab, show.zeros = FALSE, 
                              show.margins = TRUE, ... )
  {
    obj <- x
    tmp <- as.data.frame(x)
    x <- tmp[,1]
    y <- tmp[,2]
    z <- tmp[,3]
    tableflag <- TRUE

    if(missing(xlab)) xlab <- names(dimnames(obj))[1]
    if(missing(ylab)) ylab <- names(dimnames(obj))[2]
    if(missing(zlab)) zlab <- "Freq"

    balloonplot.default(x, y, z, xlab=xlab, ylab=ylab, zlab=zlab, 
                        show.zeros = show.zeros, 
                        show.margins = show.margins, ...)
  }



balloonplot.default <- function(x,y,z,
                                xlab,
                                ylab,
                                zlab=deparse(substitute(z)),
                                dotsize=2/max(strwidth(19),strheight(19)),
                                dotchar=19,
                                dotcolor="skyblue",
                                text.size=1,
                                text.color=par("fg"),
                                main,
                                label=TRUE,
                                label.digits=2,
                                label.size=1,
                                label.color=par("fg"),
                                scale.method=c("volume","diameter"),
                                scale.range=c("absolute","relative"),
                                colsrt=par("srt"),
                                rowsrt=par("srt"),
                                colmar=1,
                                rowmar=2,
                                show.zeros=FALSE,
                                show.margins=TRUE,
                                cum.margins=TRUE,
                                sorted=TRUE,
                                label.lines=TRUE,
                                fun=function(x)sum(x,na.rm=T),
                                hide.duplicates=TRUE,
                                ... )
{

  if(is.null(names(x)))
    {
      xnames <- as.character(substitute(x))
      if(length(xnames)>1) xnames <- xnames[-1]
    }
  else
     xnames <- names(x)

  if(is.null(names(y)))
    {
      ynames <- as.character(substitute(y))
      if(length(ynames)>1) ynames <- ynames[-1]
    }
  else
     ynames <- names(y)

  ####
  ## Handle arguments
  ####
  
  scale.method <- match.arg(scale.method)
  scale.range  <- match.arg(scale.range)

  if( scale.method=="absolute" && any(z < 0, na.rm=TRUE )  )
    warning("z value(s) below zero detected.",
            " No balloons will be displayed for these cells.")
  
  if(missing(main))
    {
      if(scale.method=="volume")
        main <- paste("Balloon Plot for ",
                      paste(xnames, collapse=", "),
                      " by ",
                      paste(ynames, collapse=", "),
                      ".\nArea is proportional to ", zlab, ".", sep='')
      else
        main <- paste("Balloon Plot for ",
                      paste(ynames, collapse=", "),
                      " by ",
                      paste(ynames, collapse=", "),
                      ".\nDiameter is proportional to ", zlab, ".", sep='')
      }

  if(length(dotcolor)<length(z))
    dotcolor <- rep(dotcolor, length=length(z))

  ####
  ## Make sure x and y are lists
  ####
  
  if(is.list(x))
    {
      xlabs <- x
      x$sep=":"
      x <- do.call(paste, x)
    }
  else
      xlabs <- list(x)

  if(is.list(y))
    {
      ylabs <- y
      y$sep=":"
      y <- do.call(paste, y)
      ylab <-  paste( names(y) )      
    }
  else
    ylabs <- list(y)


  ####
  ## sort everything into a useful order
  ####
  if(sorted)
    {
      ord.x <- do.call(order, xlabs)
      ord.y <- do.call(order, ylabs)
    }
  else
    ord.x <- ord.y <- 1:length(x)


  forceOrder <- function(X, sord, lord)
    factor(X[sord], levels=unique(X[lord]))
  
  x <- forceOrder(x, ord.y, ord.y)
  y <- forceOrder(y, ord.y, ord.y)
  z <- as.numeric(z[ord.y])
  dotcolor <- dotcolor[ord.y]

  xlabs <- unique(data.frame(lapply(xlabs, forceOrder,
                                    sord=ord.y, lord=ord.y)))
  ylabs <- unique(data.frame(lapply(ylabs, forceOrder,
                                    sord=ord.y, lord=ord.y)))

  ####
  ## Function to scale circles to fill the containing box
  ####
  myscale <- function(X, min=0, max=16, scale.method, scale.range)
    {
      if(scale.method=="volume")
        {
          X[X<0] <- 0
          X <- sqrt(X)
        }

      if(scale.range=="relative")
          X <- (X-min(X, na.rm=TRUE)) # put min to 0
      X <- X / max(X, na.rm=TRUE )    # put max to 1
      X <- min + X  * (max - min)     # now to [min,max]

      cin.x <- par("cin")[1]
      cin.y <- par("cin")[2]
      if(cin.x < cin.y) X <- X * cin.x/cin.y
      X
    }

  nlabels.y <- length(ylabs)
  nlabels.x <- length(xlabs)




  ####
  ## Combine duplicate entries
  ####
  # Do twice, once for data, once for colors

  tab1 <- split( data.frame(z,dotcolor,x,y), f=list(x,y) )
  ztab <- do.call(rbind,
                  lapply(
                         tab1,
                         FUN=function(X) cbind(z=fun(X[,1]),X[1,-1])
                         )
                  )
  ####
  ## Do the plotting
  ###

  oldpar <- par("xpd","mar")
  on.exit( par(oldpar) )
  #par(xpd=NA, mar=c(1,1,5,1)+0.1)   # clip drawing to device region

  if(!show.margins)
    {
      xlim=c(-0.5,nlevels(x)+nlabels.y*rowmar-0.25)   # extra space on either
                                                      # end of plot for labels
      ylim=c(0.50,nlevels(y)+nlabels.x*colmar+1)      # and so dots don't cross
                                                      # into margins,
    }
  else
    {
      xlim=c(-0.5,nlevels(x)+nlabels.y*rowmar+1)      # extra space on either
                                                      # end of plot for labels
      ylim=c(0,nlevels(y)+nlabels.x*colmar+1)         # and so dots don't cross
                                                      # into margins,
    }

  
  plot(x=nlabels.y*rowmar+0.25 + as.numeric(ztab$x) - 1,
       y=nlevels(y) - as.numeric(ztab$y) + 1,
       cex=myscale(ztab$z, max=dotsize, scale.method=scale.method, scale.range=scale.range),
       pch=dotchar, # plot character
       col=as.character(ztab$dotcolor), # dot color
       xlab="",
       ylab="",
       xaxt="n", # no x axis lables
       yaxt="n", # no y axis lables
       bty="n",  # no box around the plot
       xaxs = "i",
       yaxs = "i",
       xlim=xlim,
       ylim=ylim,
       ...
     )

  ny <- nlevels(ztab$y)
  nx <- nlevels(ztab$x)


  sumz    <- sum(ztab$z, na.rm=TRUE)
  colsumz <- sapply(split( ztab$z, ztab$y), sum, na.rm=TRUE) # works
  rowsumz <- sapply(split( ztab$z, ztab$x), sum, na.rm=TRUE) # broken
  
  if(show.margins)
    {
      ## column totals
      text(
           x=(1:nx) + nlabels.y*rowmar + 0.25 -1,
           y=0.25,
           labels=format(c(sumz, rowsumz), digits=label.digits)[-1],
           font=1,
           adj=c(0.5,0.0),
           col=text.color,
           cex=text.size
           )

      ## row totals
      rowlabs <- format(c(sumz, colsumz), digits=label.digits)[-1]
      width <- max(strwidth(rowlabs),na.rm=TRUE)
      text(
           x=nx + nlabels.y*rowmar-0.25+width,
           y= (ny:1),
           labels=rowlabs,
           font=1,
           adj=c(1.0,0.5),
           col=text.color,
           cex=text.size           
           )

      ## overall total
      text(
           x=nx + nlabels.y*rowmar-0.25+width,
           y=0.25,
           labels=sumz,
           font=1,
           adj=c(1.0,0.0),
           col=text.color,
           cex=text.size           
           )
    }
     
  if(cum.margins)
    {
      ## Row Sums at left
      cx <- c(0, cumsum(rowsumz) / sumz)
      rect(xleft   = nlabels.y*rowmar - 1 - 0.25 + 1:nx,
           xright  = nlabels.y*rowmar - 1 + 0.75 + 1:nx,
           ybottom = ny+0.75+cx[1:nx]*colmar*nlabels.x,
           ytop    = ny+0.75+cx[2:(nx+1)]*colmar*nlabels.x,
           col     = "lightgray",
           border  = NA)

      ## Col Sums at top
      cy <- c(0, cumsum(colsumz) / sumz)
      rect(xleft   = -0.5 +rowmar*cy[ny:1]*nlabels.y,
           xright  = -0.5 +rowmar*cy[(ny+1):2]*nlabels.y,
           ybottom = 1:ny-0.5,
           ytop    = 1:ny+0.5,
           col     = "lightgray",
           border  = NA)
      
      
      tx <- paste(levels(x),"\n[",rowsumz,"]")
      ty <- paste(levels(y),"\n[",colsumz,"]")
    }

  
  ###
  ## Horizontal borders between cells
  ###
  segments(
           x0=nlabels.y*rowmar-0.25,
           x1=nx+nlabels.y*rowmar-0.25,
           y0=(0:ny)+0.5,
           y1=(0:ny)+0.5
           )
  
  ###
  ## Vertical borders between cells
  ###
  segments(
           x0=(0:nx)+nlabels.y*rowmar-0.25,
           x1=(0:nx)+nlabels.y*rowmar-0.25,
           y0= 0.5,
           y1=ny+0.5,
           )


  if(hide.duplicates)
    undupe <- function(X) 
      {
                                        # convert duplicates into blanks
        X <- as.character(X)
        c(X[1], ifelse(X[-1] == X[-length(X)], "", X[-1]))
      }
  else
    undupe <- function(X) X

  ### 
  ## Column labels
  ###
  for(i in 1:nlabels.x)
    {
      y <- ny + 0.75 + (nlabels.x - i + .5)*colmar
      text(
           x= (1:nx) + nlabels.y*rowmar + 0.25 - 1,
           y= y,
           labels=undupe(xlabs[,i]),
           srt=colsrt,
           font=1,
           col=text.color,
           cex=text.size
           )
    }

  ### 
  ## Row labels
  ###
  for(i in 1:length(ylabs))
    {
      text(
           y=ny:1,
           x= (i-0.5)*rowmar-0.5,
           labels=undupe(ylabs[,i]),
           srt=rowsrt,
           font=1,
           col=text.color,
           cex=text.size
           )
    }

  ####
  ## Column headers for row labels
  ####
  if(missing(ylab) || length(ylab)==0)
    text(
         x=((1:length(ylabs))-0.5)*rowmar-0.5,
         y=ny+0.5,
         labels=ynames,
         srt=colsrt,
         font=2,
         adj=c(0.5,0.0),
         col=text.color,
         cex=text.size
         )
  else
    text( 
         x=((1:length(ylab))-0.5)*rowmar-0.5,
         y=ny+0.5,
         labels=ylab,
         srt=colsrt,
         font=2,
         adj=c(0.5,0.0),
         col=text.color,
         cex=text.size
         )

  ####
  ## Row headers for column labels
  ####
  if(missing(xlab) || length(xlab)==0)
    text(
         x= nlabels.y*rowmar - 0.25 - strwidth(','),
         y= ny + 0.75 + ((nlabels.x:1) - 1 + .5)*colmar,
         labels=xnames,
         srt=colsrt,
         font=2,
         adj=c(1,0.5),
         col=text.color,
         cex=text.size
         )
  else
    text(
         x= nlabels.y*rowmar - 0.25 - strwidth(','),
         y= ny + 0.75 + ((length(xlab):1) - 1 + .5)*colmar,
         labels=xlab,
         srt=colsrt,
         font=2,
         adj=c(1,0.5),
         col=text.color,
         cex=text.size
         )

  ###
  ## add borders to row and column headers
  ###
  if(label.lines)
    {
      segments(                          # left: vertical lines
               x0=(0:nlabels.y)*rowmar-0.5,
               x1=(0:nlabels.y)*rowmar-0.5,
               y0=0.5,
               y1=ny+0.5
               )
      
      segments(
               x0=nlabels.y*rowmar-0.25,        # top: horizontal lines
               x1=nlabels.y*rowmar + nx - 0.25,
               y0=(0:nlabels.x)*colmar  +ny+0.75,
               y1=(0:nlabels.x)*colmar  +ny+0.75
               )
    }


  ####
  ## annotate cells with actual values
  ####
  if(label){
    if(show.zeros) 
     indiv <- 1:length(ztab$y) 
    else 
      indiv <- which(ztab$z != 0)
    
    text(x=as.numeric(ztab$x[indiv])+ nlabels.y*rowmar - 0.75,     # as.numeric give numeric values
         y=ny - as.numeric(ztab$y[indiv]) + 1,
         labels=format(ztab$z[indiv], digits=label.digits),       # label value
         font=2,
         adj=c(0.5,0.5),
         col=label.color,
         cex=label.size
         )
  }
  # put a nice title
  title(main=main)
}
