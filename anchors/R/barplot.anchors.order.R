## Plot frequency of vignette ordering
## 
## Author:    Jonathan Wand <wand(at)stanford.edu>
## Created:   2006-10-02
## Modified:  $Date: $
## Revision:  $Revision:  $
## RCS-ID:    $Id: $

barplot.anchors.order <- function(height,...,top = 20) {

  op <- par(no.readonly=TRUE)
  if (height$ties == "interval") {
    par(mar= par()$mar + c(0,1,0,0) )
  }
  
  args <- list(...)

  
  fn <- function(x,y) {
    if (!is.null(x))
      return(x)
    else
      return(y)
  }
  
  mf <- match.call()
  mf$x <- mf$top <- NULL
  mf[[1]] <- as.name("barplot.default")

  mf$cex.names <- fn( args$cex.names, .8)
  mf$horiz     <- fn( args$horiz    , TRUE)
  mf$las       <- fn( args$las      , 1   )
  mf$axes      <- fn( args$axes     , TRUE)
  mf$main      <- fn( args$main     , height$main )
  mf$height    <- height$freq[1:top]

  mf$ylim      <- fn( args$ylim     , c(0, min(top,length(height$freq)) ) )
  mf$space     <- fn( args$space    , .2)
  mf$width     <- fn( args$width    , .8)

#  mf$ylim      <- fn( args$ylim     , c(0,1) )
#  mf$width     <- 1/length(height$freq)
  
  
  mf$xlab      <- fn( args$xlab     , "Frequency" )
  
  bp <- eval(mf)
  par(op)
  
  return(invisible(bp))

}
