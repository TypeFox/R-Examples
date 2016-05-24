#' Perform Principal Component Analysis on a large data set
#'
#' Run \code{prcomp} on subsamples of the data set and compile the results for the first dimension.
#' 
#' @param x data.frame, data over which to run PCA
#' @param filename character, name of the file containing the data. This must be a tab-delimited file with a header row formatted per the default options for \code{\link{read.delim}}.
#' @param db Object type, database connection to table containing the data (NOT IMPLEMENTED).
#' @param subsample numeric or logical, If an integer, size of each subsample.  If FALSE, runs PCA on entire data set.
#' @param n.subsamples numeric, number of subsamples.
#' @param ignore.cols numeric, indices of columns not to include.
#' @param use.cols numeric, indices of columns to use.
#' @param return.sds logical, if TRUE return the standard deviations of each network's edge weights.
#' @param progress.bar logical, if TRUE then progress in running subsamples will be shown.
#' @return If \code{return.sds} is FALSE, return named vector of 
#' component weights for first dimension of principal component 
#' analysis (see example for comparison to \code{\link{prcomp}}).
#' 
#' If \code{return.sds} is TRUE, return a list.
#' \tabular{ll}{
#' coefficients \tab named vector of the component weights for first dimension of principal component analysis (see example for comparison to \code{\link{prcomp}}).\cr
#' sds \tab named vector of the standard deviations of each network's edge weights.\cr
#' }
#' @details Scales the function \code{\link{prcomp}} to data sets with an arbitrarily
#' large number of rows by running \code{prcomp} on repeated subsamples of the rows. 
#' @export
#' @seealso \code{\link{prcomp}}
#' @references
#' \url{https://github.com/shaptonstahl/}
#' @author Stephen R. Haptonstahl \email{srh@@haptonstahl.org}
#' @examples
#' data(iris)        # provides example data
#' prcomp(iris[,1:4], center=FALSE, scale.=FALSE)$rotation[,1]
#' ScalablePCA(iris, subsample=10, use.cols=1:4)
#' ScalablePCA(iris, subsample=10, ignore.cols=5)
ScalablePCA <- function(x,
                        filename=NULL,
                        db=NULL,
                        subsample=1e4,
                        n.subsamples=1e3,
                        ignore.cols,
                        use.cols,
                        return.sds=FALSE,
                        progress.bar=FALSE) {
  # Guardians
  if( missing(x) ) {
    if( is.null(filename) ) {
      if( is.null(db) ) {
        stop("One of x, filename, or db must be specified.")
      } else {
        source.type <- "db"
        # ADD CODE: TEST CONNECTION TO DATABASE
        stop("db option not yet implemented")
      }
    } else {
      source.type <- "file"
      if(length(filename) > 1) {
        filename <- filename[1]
        warning("filename has multiple values; using only the first value")
      }
      if( is.character(filename) ) {
        if( !file.exists(filename) ) {
          stop("Unable to make a connection to ", filename)
        }
      } else {
        stop("filename must be a character (length=1) vector specifying the name of the file to open")
      }
    }
  } else {
    source.type <- "data.frame"
    if( !is.data.frame(x) ) {
      stop("x must be a data.frame")
    }
  }
  
  subsample.error.message <- "subsample must be a positive whole number or FALSE"
  if(length(subsample) > 1) {
    subsample <- subsample[1]
    warning("subsample has multiple values; using only the first value")
  }
  if( is.logical(subsample) ) {
    if( subsample ) {
      stop(subsample.error.message)
    } else {
      # perform PCA on entire data set
      n.subsamples <- 1
    }
  } else {
    if( is.numeric(subsample) ) {
      if( subsample <= 0 ) stop(subsample.error.message)
      if( 0 != subsample %% 1 ) stop(subsample.error.message)
      
    } else {
      stop(subsample.error.message)
    }
  }

  n.subsamples.error.message <- "n.subsamples must be a positive whole number"
  if(length(n.subsamples) > 1) {
    n.subsamples <- n.subsamples[1]
    warning("n.subsamples has multiple values; using only the first value")
  }
  if( is.numeric(n.subsamples) ) {
    if( n.subsamples <= 0 ) {
      stop(n.subsamples.error.message)
    }
  } else {
    stop(n.subsamples.error.message)
  }
  
  ignore.cols.error.message <- "ignore.cols must be numeric with positive whole values"
  if( !missing(ignore.cols) ) {
    ignore.cols <- sort(unique(ignore.cols))
    if( !is.numeric(ignore.cols) ) stop(ignore.cols.error.message)
    if( any(ignore.cols <= 0) ) stop(ignore.cols.error.message)
    if( any(0 != ignore.cols %% 1) ) stop(ignore.cols.error.message)
  }
  use.cols.error.message <- "use.cols must be numeric with positive whole values"
  if( !missing(use.cols) ) {
    use.cols <- sort(unique(use.cols))
    if( !is.numeric(use.cols) ) stop(use.cols.error.message)
    if( any(use.cols <= 0) ) stop(use.cols.error.message)
    if( any(0 != use.cols %% 1) ) stop(use.cols.error.message)
  }
  # End Guardians
  
  # Set internal sampling function to one of three externally-defined functions
  if( "db" == source.type ) {
    # n.rows.full.dataset <- ADD CODE TO READ db AND GET NUMBER OF RECORDS
    pca.sample <- function(n=subsample) GetSampleFromDb(n=n, db=db)
  } else if( "file" == source.type ) {
    # Count rows in the text file, -1 for the header
    n.rows.full.dataset <- -1  # don't count header row
    con <- file(filename, "r")
    while( length(input <- readLines(con, n=1e4)) > 0 ) {
      n.rows.full.dataset <- n.rows.full.dataset + length(input)
    }
    close(con)
    
    pca.sample <- function(n=subsample) GetSampleFromFile(n=n, 
                                                          out.of=n.rows.full.dataset, 
                                                          filename=filename)
  } else if( "data.frame" == source.type ) {
    n.rows.full.dataset <- nrow(x)
    pca.sample <- function(n=subsample) GetSampleFromDataFrame(n=n, x=x)
  } else{
    stop("source.type must be one of db, file, or data.frame")
  }
  
  if( missing(x) && identical(subsample, FALSE) ) subsample <- n.rows.full.dataset
  
  # prepare matrix to store returned weights
  samp <- pca.sample(n=1)
  if( !missing(ignore.cols) ) {
    use.cols <- setdiff(1:ncol(samp), ignore.cols)
  } else {
    if( !all(use.cols %in% 1:ncol(samp)) ) {
      stop("use.cols must set values between 1 and ", ncol(samp), " for this data set")
    }
  }
  n.cols.used <- length(use.cols)
  
  if( "db" == source.type ) {
    # ADD CODE: SET use.names TO NAMES OF COLUMNS SUBSETTED USING use.cols
  } else if( "file" == source.type ) {
    # ADD CODE: SET use.names TO NAMES OF COLUMNS SUBSETTED USING use.cols
  } else if( "data.frame" == source.type ) {
    use.names <- names(x)[use.cols]
  } else{
    stop("source.type must be one of db, file, or data.frame")
  }
  
  draws <- matrix(0, nrow=n.subsamples, ncol=n.cols.used)
  colnames(draws) <- use.names
  
  if(return.sds) {
    sd.draws <- matrix(0, nrow=n.subsamples, ncol=n.cols.used)
    colnames(sd.draws) <- use.names
  }
  
  # perform the function
  if(progress.bar) pb <- txtProgressBar(max=n.subsamples, style=2)
  for(g in 1:n.subsamples) {
    samp <- pca.sample()
    pca.result <- prcomp(samp[,use.cols], center=FALSE, scale.=FALSE)
    draws[g,] <- pca.result$rotation[,1]
    
    if(return.sds) {
      sd.draws[g,] <- apply(samp[,use.cols], 2, sd)
    }
    
    if(progress.bar) setTxtProgressBar(pb, g)
  }
  if(progress.bar) close(pb)
    
  # prepare and return the output
  out <- colMeans(draws)
  out <- abs(out)
  if(return.sds) {
    out <- list(coefficients=out,
                sds=colMeans(sd.draws))
  }
  return(out)
}