#' This is a utility function to compare two covariance matrices
#' 
#' @details
#' This method takes in two different covariance/correlation matrices computed
#' using two different methods and visullay compares them using the ellipse plot.
#' It produces a matrix with ellipses drawn in the upper triangle. The ellipse 
#' is drawn to be a contour of a standard bivariate normal with correlation given
#' by the correlation of the two assets. One ellipse is drawn for each covariance
#' matrix.
#' 
#' @import robust
#' @param  cov1 covariance matrix using the first method
#' @param  cov2 covariance matrix using the second method
#' @param  labels strings indicating the type of methods used for comparison
#' @param  corr flag indicating if the supplied matrices are of type covariance
#'          or correlation
#' @author Rohit Arora
#' @export
#' 
compareCov <- function(cov1, cov2, labels, corr=FALSE) {
  
    if (length(labels) != 2) stop("There must be two labels")
    
    if (!all(dim(cov1) == dim(cov2))) stop("Matrix dimensions are unequal")
    if(nrow(cov1) != ncol(cov1)) stop("Matrix is not square")
    if(nrow(cov1) == 2) stop("Need more than 2 dims")
  
    complist <- list(list(corr=corr, cov=cov1), 
                     list(corr=corr, cov=cov2))
    names(complist) <- labels
    ellipsesPlot.covfm(complist)
}

#' Plot data to visualize missing values
#' 
#' @details
#' This method takes in data as an xts object and plots the data. 
#' Missing values highlighted in red for matrix plot and time series of returns 
#' are shown in in Summary plot
#' 
#' @param  data an xts/zoo object
#' @param  which takes values 3/4. 3 = Summary plot, 4 = Matrix plot
#' @import VIM reshape2 ggplot2 xts zoo
#' @author Rohit Arora
#' @export
#' 
#' 
plotmissing <- function(data, which=c(3,4)) {
  cols <- colnames(data)
  if (length(cols) == 0) stop("Data should have column names")
  
  which <- which[1]
  
  options(warn=-1)
  
  if (which == 3)  {
    ind <- which.min(apply(is.na(data),  2, which.min))
    dates <- index(data[,ind])
    
    d <- melt(coredata(data)); colnames(d) <- c("Index","Symbol","Returns")
    symCount <- ncol(data)
    
    year.dates <- format(dates, "%Y")
    ind <- sapply(unique(year.dates), function(val) 
      which.max(year.dates == val))
    ind.ind <- seq.int(1, length(ind), length.out = min(10,length(ind)))
    ind <- ind[ind.ind]
    
    p <- ggplot(data=d, aes_string(x='Index', y='Returns', colour='Symbol', group='Symbol')) + 
      geom_line() + xlab("Dates") + ylab("Returns") +
      scale_x_discrete(breaks = ind, labels=year.dates[ind]) +
      facet_wrap(~Symbol, ncol=round(sqrt(symCount)), scales = "free_x") + 
      theme(legend.position="none")
    
    print(p)
  }
  
  if(which == 4)  {
    data <- coredata(data)
    
    if(class(data) == "data.frame") cols <- cols[unlist(lapply(data, is.numeric))]
    data <- data[, cols]
    colnames(data) <- sapply(cols, function(name) substr(name,1,9))
    matrixplot(data, main="Location of Missing values")
  }
  
  options(warn=0)
}