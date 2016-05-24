# $Id: plotCI.R 1948 2015-04-23 21:23:36Z warnes $


plotCI <- function (x,
                    y = NULL,
                    uiw,
                    liw = uiw,
                    ui,
                    li,

                    err='y',
                    ylim=NULL,
                    xlim=NULL,
                    type="p",

                    col=par("col"),
                    barcol=col,
                    pt.bg = par("bg"),

                    sfrac = 0.01,
                    gap=1,

                    lwd=par("lwd"),
                    lty=par("lty"),

                    labels=FALSE,

                    add=FALSE,
                    xlab,
                    ylab,

                    minbar,
                    maxbar,
                    ...
                    )
{

  if (is.list(x)) {
    y <- x$y
    x <- x$x
  }

  if(invalid(xlab))
    xlab <- deparse(substitute(x))

  if(invalid(ylab))
    {
      if(is.null(y))
        {
          xlab  <- ""
          ylab <- deparse(substitute(x))
        }
      else
        ylab <- deparse(substitute(y))
    }

  if (is.null(y)) {
    if (is.null(x))
      stop("both x and y NULL")
    y <- as.numeric(x)
    x <- seq(along = x)
  }


  if(err=="y")
    z  <- y
  else
    z  <- x

  if(invalid(uiw))
    uiw <- NA
  if(invalid(liw))
    liw <- NA
  
  
  if(invalid(ui))
    ui <- z + uiw
  if(invalid(li))
    li <- z - liw

  if(!invalid(minbar))
    li <- ifelse( li < minbar, minbar, li)

  if(!invalid(maxbar))
    ui <- ifelse( ui > maxbar, maxbar, ui)

   if(err=="y")
     {
       if(is.null(ylim))
         ylim <- range(c(y, ui, li), na.rm=TRUE)
       if(is.null(xlim) && !is.R() )
         xlim <- range( x, na.rm=TRUE)
     }
   else if(err=="x")
     {
       if(is.null(xlim))
         xlim <- range(c(x, ui, li), na.rm=TRUE)
       if(is.null(ylim) && !is.R() )
         ylim <- range( x, na.rm=TRUE)
     }

  if(!add)
    {
      if(invalid(labels) || labels==FALSE )
        plot(x, y, ylim = ylim, xlim=xlim, col=col,
             xlab=xlab, ylab=ylab, type="n", ...)
      else
        {
          plot(x, y, ylim = ylim, xlim=xlim, col=col, type="n",
               xlab=xlab, ylab=ylab,  ...)
          text(x, y, label=labels, col=col, ... )
        }
    }
  if(err=="y")
    {
      if(gap!=FALSE)
        gap <- strheight("O") * gap
      smidge <- par("fin")[1] * sfrac


      # draw upper bar
      if(!is.null(li))
          arrows(x , li, x, pmax(y-gap,li), col=barcol, lwd=lwd,
                 lty=lty, angle=90, length=smidge, code=1)
      # draw lower bar
      if(!is.null(ui))
          arrows(x , ui, x, pmin(y+gap,ui), col=barcol,
                 lwd=lwd, lty=lty, angle=90, length=smidge, code=1)
    }
  else
    {
      if(gap!=FALSE)
        gap <- strwidth("O") * gap
      smidge <- par("fin")[2] * sfrac

      # draw left bar
      if(!is.null(li))
        arrows(li, y, pmax(x-gap,li), y, col=barcol, lwd=lwd,
                 lty=lty, angle=90, length=smidge, code=1)
      if(!is.null(ui))
        arrows(ui, y, pmin(x+gap,ui), y, col=barcol, lwd=lwd,
                 lty=lty, angle=90, length=smidge, code=1)

    }

  ## _now_ draw the points (to avoid having lines drawn 'through' points)
  points(x, y, col = col, lwd = lwd, bg = pt.bg, type = type, ...)

  invisible(list(x = x, y = y))
}
