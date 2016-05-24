# $Id: boxplot2.R 2052 2015-06-02 19:03:20Z warnes $
boxplot.n <- function( ..., top=FALSE, shrink=1.0, textcolor=NULL )
    {
        .Defunct("boxplot2", package="gplots")
    }

boxplot2  <- function( ..., top=FALSE, shrink=1.0, textcolor=NULL )
  {
    boxcall <- match.call()           # get call
    boxcall$top  <- boxcall$shrink  <- boxcall$textcolor  <- NULL
    boxcall[[1]]  <- as.name("boxplot")

    if(is.R())
      {
        box <- eval(boxcall, parent.frame())
        mids <- 1:length(box$n)
      }
    else
      {
        mids <- eval(boxcall, parent.frame())
        boxcall$plot <- FALSE
        box <- eval(boxcall, parent.frame())
      }

    if(top)
      {
        where  <- par("usr")[4]
        adj  <- c(0.5,1)
      }
    else
      {
        where  <- par("usr")[3]
        adj  <- c(0.5,0)
      }
    tcex <- par("cex")
    par(cex=shrink*tcex)

    if(is.R())
      text( x=mids, y=where, labels=paste("n=",box$n,sep=""), adj=adj,
			      col=textcolor)
    else
      {
        if( is.null(textcolor) )
          textcolor <- 1
        space <- ifelse(top, -1, 1) * par("1em")[2] / 2

        text( x=mids, y=where + space, labels=paste("n=",box$n,sep=""), adj=adj[1],
             col=textcolor)

      }

    par(cex=tcex)

    invisible(box)
  }

