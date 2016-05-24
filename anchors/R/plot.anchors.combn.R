#######################################################################
##
## Function: plot.anchors.combn()
## Author  : Jonathan Wand <wand@stanford.edu>
## 
#######################################################################
plot.anchors.combn <- function(x, ...,
                               xy = c("minimum","interval")) {

  if (!any( "anchors.combn" == class(x)) )
    stop("x must be of class 'anchors.combn'")

  idx <- pmatch( tolower(xy), tolower(colnames(x)), 0)
  
  if (length(idx) != 2 || !all(idx > 0) )
    stop("'xy' argument must specify names of two (and only two) columns of the anchors.combn object")

  pts   <- na.omit(x)
  args  <- list(...)
  nsets <- length(table(nchar(pts$Vignettes)))

  cat("Requested plot:\n")
  cat(" x-axis: ", colnames(pts)[idx[1]], "\n")
  cat(" y-axis: ", colnames(pts)[idx[2]], "\n")
#  cat(" using indices", idx, "\n")

  
  ## coloring of sets of vignettes
  if (!is.null(args$col)) {
    if (length(args$col) !=  nsets) {
      col <- rep(args$col,nsets)[1:nsets]
    } else {
      col <- args$col
    }
  } else {
    col <- rep( c("black","grey"), nsets)[1:nsets]
  }

    fn <- function(x,y) {
      if (!is.null(x))
        return(x)
      else
        return(y)
    }

  
  mf <- match.call()
  m  <- match( c("xlim","ylim","log","main","sub","xlab","ylab","ann","axes",
                 "frame.plot","panel.first","panel.last","asp"), names(mf),0)
  mf <- mf[ c(1,m) ]
  
  mf[[1]] <- as.name("plot.default")

  mf$xlab       <- fn( args$xlab , colnames(pts)[ idx[1] ] )
  mf$ylab       <- fn( args$xlab , colnames(pts)[ idx[2] ] )
  mf$xlim       <- fn( args$xlim , range(pts[ ,idx[1]] ))
  mf$ylim       <- fn( args$ylim , range(pts[ ,idx[2]] ))
  mf$frame.plot <- fn( args$frame.plot, FALSE )
  mf$x <- 0
  mf$y <- 0
  mf$type <- "n"
  
  ## do plot()
  eval(mf)

  ## do text()
  mft <- match.call()
  m   <- match( c("adj","pos","vfont","cex","font"), names(mft), 0)
  mft <- mft[ c(1,m) ]
  mft[[1]] <- as.name("text.default")

  for (i in 1:nrow(pts)) {
    txt <- as.character(pts$Vignettes[i])
    mft$labels <- txt
    mft$x      <- pts[i,idx[1]]
    mft$y      <- pts[i,idx[2]]
#    cat("observation",mft$labels, mft$x, mft$y, "\n")
    mft$col    <- col[nchar(txt)]
    eval(mft)
  }
  
  return(invisible(pts[,c(1,idx)]))
}

