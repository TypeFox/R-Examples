svm.data.root <- system.file("data_sets", "svm", package="gmum.r")
svm.colon_cancer.path <- file.path(svm.data.root, "colon-cancer")

svm.lib.libsvm <- "libsvm"
svm.lib.svmlight <- "svmlight"

svm.prep.2e <- "2e"
svm.prep.none <- "none"

svm.kernel.linear <- "linear"
svm.kernel.poly <- "poly"
svm.kernel.rbf <- "rbf"
svm.kernel.sigmoid <- "sigmoid"

svm.plot.contour <- "contour"
svm.plot.pca <- "pca"

read.libsvm = function( filename, dimensionality ) {
  
  content = readLines( filename )
  num_lines = length( content )
  yx = matrix( 0, num_lines, dimensionality + 1 )
  
  # loop over lines
  for ( i in 1:num_lines ) {
    
    # split by spaces
    line = as.vector( strsplit( content[i], ' ' )[[1]])
    
    # save label
    yx[i,1] = as.numeric( line[[1]] )
    
    # loop over values
    for ( j in 2:length( line )) {
      
      # split by colon
      index_value = strsplit( line[j], ':' )[[1]]
      
      index = as.numeric( index_value[1] ) + 1  	# +1 because label goes first
      value = as.numeric( index_value[2] )
      
      yx[i, index] = value
    }
  }
  
  return( yx )
}

svm.dataset.colon_cancer <- function() {
  bc <- read.libsvm(svm.colon_cancer.path, 2000)
  return(bc)
}

svm.dataset.circles <- function() { 
  matrix( 
    c(0,1,0,1,0,0,1,1,0,1,1,0),
    ncol=3,
    nrow=4,
    dimnames=list(c(),c("x","y","t")))
}

#' @title Measure accuracy scoreof a prediction
#' 
#' @description Calculates accuracy of a prediction, returns precent of correctly predicted examples 
#' over all test examples.
#' @export svm.accuracy
#' @rdname svm.accuracy
#' 
#' @usage svm.accuracy(prediction, target)
#' 
#' @param prediction factor or 1 dim vector with predicted classes
#' @param target  factor or 1 dim vector with true classes
#' 
#' #' @examples
#' \dontrun{
#' # firstly, SVM model needs to be trained
#' svm <- SVM(x, y, core="libsvm", kernel="linear", C=1)
#' # then we can use it to predict unknown samples
#' p <- predcit(svm, x_test)
#' acc <- svm.accuracy(p, y)
#' }
svm.accuracy <- function(prediction, target) {
    if ( length(target) != length(prediction)) {
      stop("Prediction's and target's length don't match!")
    }

    diff = as.numeric(as.factor(target)) -  as.numeric(as.factor(prediction))
    acc <- sum(diff == 0) / length(target)
    return(acc) 
}
 
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

scale.data.frame <-
  function(x, center = TRUE, scale = TRUE)
  {
    i <- sapply(x, is.numeric)
    if (ncol(x[, i, drop = FALSE])) {
      x[, i] <- tmp <- scale.default(x[, i, drop = FALSE], na.omit(center), na.omit(scale))
      if(center || !is.logical(center))
        attr(x, "scaled:center")[i] <- attr(tmp, "scaled:center")
      if(scale || !is.logical(scale))
        attr(x, "scaled:scale")[i]  <- attr(tmp, "scaled:scale")
    }
    x
  }

read.matrix.csr <- function(file, fac = TRUE, ncol = NULL) {
  l <- strsplit(readLines(file), "[ ]+")
  
  ## extract y-values, if any
  y <- if (is.na(l[[1]][1]) || length(grep(":",l[[1]][1])))
    NULL
  else
    sapply(l, function(x) x[1])
  
  ## x-values
  rja <- do.call("rbind",
                 lapply(l, function(x)
                   do.call("rbind",
                           strsplit(if (is.null(y)) x else x[-1], ":")
                   )
                 )
  )
  ja <- as.integer(rja[,1])
  ia <- cumsum(c(1, sapply(l, length) - !is.null(y)))
  
  max.ja <- max(ja)
  dimension <- c(length(l), if (is.null(ncol)) max.ja else max(ncol, max.ja))
  x = new(getClass("matrix.csr", where = asNamespace("SparseM")),
          ra = as.numeric(rja[,2]), ja = ja,
          ia = as.integer(ia), dimension = as.integer(dimension))
  if (length(y))
    list(x = x, y = if (fac) as.factor(y) else as.numeric(y))
  else x
}
