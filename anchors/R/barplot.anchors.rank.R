#######################################################################
##
## Function: barplot.anchors()
## Author  : Jonathan Wand <wand@stanford.edu>
##
#######################################################################
barplot.anchors.rank <- function(height, ... ,
                                 ties = c("uniform","minentropy","omit","cpolr") ) {

  debug <- 0
  
  ties <- match.arg(ties)

  if ((class(height) != "anchors.rank"))
    stop("anchors.rank objects must be listed first in function")
  
  args <- list(...)

    fn <- function(x,y) {
      if (!is.null(x))
        return(x)
      else
        return(y)
    }


  unconditional <- fn( args$unconditional,  FALSE)
  
  mf <- match.call()
  mf[[1]] <- as.name("barplot.default")
  mf$debug <- mf$ties <- mf$height <- mf$x <- mf$y <- mf$unconditional <- NULL

  mf$beside        <- fn( args$beside   , TRUE) 
  mf$main          <- fn( args$main     , ties )                     
  mf$names.arg     <- fn( args$names.arg, 1:height$rank$max)       
  mf$xlab          <- fn( args$xlab     , height$type)                  
  mf$ylab          <- fn( args$ylab     , "Proportion")            

  
  w <- NULL
  w <-  rbind(w,
         fitted.anchors.rank( height, ties=ties, average=TRUE, unconditional=unconditional))
  if (length(args) > 0) {
    for (i in 1:length(args)) {

      if ((class(args[[i]]) != "anchors.rank"))
        break;
      
      w <- rbind(w,
                 fitted.anchors.rank( args[[i]], ties=ties, average=TRUE, unconditional=unconditional))
      ## keep popping the second item off the stack
      mf[[2]] <- NULL
    }
  }

#  mf$xlim      <- fn( args$xlim     , c(0, height$rank$max ) )
#  mf$space     <- fn( args$space    , .1)
#  mf$width     <- fn( args$width    , .9 / nrow(w) )
  
  mf$height <- w
  
  if (debug > 0) {
    print(w)
    cat("do eval\n")
    print(mf)
  }
  
  bp <- eval(mf)

  return(invisible(bp))
}
