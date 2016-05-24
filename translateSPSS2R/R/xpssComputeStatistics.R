


#' Computes the coefficient of variation
#'
#' Helper Function for xpssCompute. R Implementation of the SPSS \code{CFVAR} Function. 
#' 
#' @usage computeCfvar(x,...)
#'
#' @param x atomic numeric or numeric vector or numeric matrix
#' @param ... further arguments passed to or from other methods.
#'
#' @return Numeric. Returns the coeffiecent of variation. 
#' If the data contains missing values, it is possible to specify an na remove command, this is \code{na.rm}. The default for na.rm is \code{na.rm=F}, if the value get changed to \code{na.rm=T} every existing missing value get omitted.
#' @author Bastian Wiessner
#' @seealso \code{\link{max}}
#' @keywords internal
#' @examples
#' xpssCompute(x = fromXPSS, variables = c("V5","V7_1"),fun = "computeCfvar")
#' 
#' @export
#' 

computeCfvar <- function(x,...){
  options(warn=-1)
  if(is.matrix(x)){
    out <- x
    val1 <- x
    val2 <- x
    for(i in 1:nrow(x)){
      val1[i,] <- mean(x[i,],...)  
      val2[i,] <- sd(x[i,],...)  
      }
      val1 <-  as.numeric(val1[,1])
      val2 <-  as.numeric(val2[,1])
    out <- val2/val1
    pos <- which(is.nan(out))
    out[pos] <- NA    
    }
  
  return(out)
}





#' Computes the Maxima
#'
#' Helper Function for xpssCompute. R Implementation of the SPSS \code{MAX} Function. 
#' 
#' @usage computeMax (x,...)
#'
#' @param x atomic numeric or numeric vector or numeric matrix
#' @param ... further arguments passed to or from other methods.
#'
#' @return Numeric. Returns the maximum value of its arguments that have valid values. 
#' If the data contains missing values, it is possible to specify an na remove command, this is \code{na.rm}. The default for na.rm is \code{na.rm=F}, if the value get changed to \code{na.rm=T} every existing missing value get omitted.
#' @author Bastian Wiessner
#' @seealso \code{\link{max}}
#' @keywords internal
#' @examples
#' xpssCompute(x = fromXPSS, variables = c("V5","V7_1"),fun = "computeMax")
#' 
#' @export
#' 

computeMax <- function(x,...){
  options(warn=-1)
  if(is.matrix(x)){
    out <- x    
    for(i in 1:nrow(x)){
      out[i,] <- max(x[i,],...)  
    }
    out <-  as.numeric(out[,1])  
    pos <- which(is.infinite(out))
    out[pos] <- NA
  }    
  options(warn=0)
  return(out)
}





#' Atithmetic Mean
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{MEAN} Function. 
#'
#' @usage computeMean (x,...)
#'
#' @param x Numeric or numeric vector
#' @param ... further arguments passed to or from other methods.
#' @return Numeric. Returns the arithmetic mean of its arguments that have valid, nonmissing values. 
#' If the data contains missing values, it is possible to specify an na remove command, this is \code{na.rm}. The default for na.rm is \code{na.rm=F}, if the value get changed to \code{na.rm=T} every existing missing value get omitted.
#' @author Bastian Wiessner
#' @seealso \code{\link{mean}}
#' @keywords internal
#' @examples
#' xpssCompute(x = fromXPSS, variables = c("V5","V7_1"),fun = "computeMean")
#' 
#' @export
#' 

computeMean <- function(x,...){
  
  if(is.matrix(x)){
    out <- x
    for(i in 1:nrow(x)){
      out[i,] <- mean(x[i,],...)  
    }
    out <-  as.numeric(out[,1])
    pos <- which(is.nan(out))
    out[pos] <- NA
  }
  return(out)
}





#' Median Value
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{MEDIAN} Function. 
#'
#' @usage computeMedian (x,...)
#'
#' @param x Numeric or numeric vector
#' @param ... further arguments passed to or from other methods.
#'
#' @return Numeric. Returns the median (50th percentile) of its arguments that have valid, nonmissing values.  
#' If the data contains missing values, it is possible to specify an na remove command, this is \code{na.rm}. The default for na.rm is \code{na.rm=F}, if the value get changed to \code{na.rm=T} every existing missing value get omitted.
#' @author Bastian Wiessner
#' @seealso \code{\link{median}}
#' @keywords internal
#' @examples
#' xpssCompute(x = fromXPSS, variables = c("V5","V7_1"),fun = "computeMedian")
#' @export

computeMedian <- function(x,...){
  if(is.matrix(x)){
    out <- x
    for(i in 1:nrow(x)){
      out[i,] <- median(x[i,],...)  
    }
    out <-  as.numeric(out[,1])
  }
  return(out)
}





#' Minima
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{MIN} Function. 
#'
#' @usage computeMin (x,...)
#'
#' @param x Numeric or numeric vector
#' @param ... further arguments passed to or from other methods.
#'
#' @return Numeric or string. Returns the minimum value of its arguments that have valid, nonmissing values. 
#' If the data contains missing values, it is possible to specify an na remove command, this is \code{na.rm}. The default for na.rm is \code{na.rm=F}, if the value get changed to \code{na.rm=T} every existing missing value get omitted.
#' @author Bastian Wiessner
#' @seealso \code{\link{min}}
#' @keywords internal
#' @examples
#' xpssCompute(x = fromXPSS, variables = c("V5","V7_1"),fun = "computeMin")
#' @export

computeMin <- function(x,...){
  options(warn=-1)
  if(is.matrix(x)){
    out <- x
    for(i in 1:nrow(x)){
      out[i,] <- min(x[i,],...)  
    }
    out <-  as.numeric(out[,1])
    pos <- which(is.infinite(out))
    out[pos] <- NA
  }
  options(warn=0)
  return(out)
}






#' Standard Deviation
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{SD} Function. 
#'
#' @usage computeSd (x,...)
#'
#' @param x Numeric or numeric vector
#' @param ... further arguments passed to or from other methods.
#'
#' @return Numeric. Returns the standard deviation of its arguments that have valid, nonmissing values.  
#' If the data contains missing values, it is possible to specify an na remove command, this is \code{na.rm}. The default for na.rm is \code{na.rm=F}, if the value get changed to \code{na.rm=T} every existing missing value get omitted.
#' @author Bastian Wiessner
#' @seealso \code{\link{sd}}
#' @keywords internal
#' @examples
#' xpssCompute(x = fromXPSS, variables = c("V5","V7_1"),fun = "computeSd")
#' @export

computeSd <- function(x,...){
  if(is.matrix(x)){
    out <- x
    for(i in 1:nrow(x)){
      out[i,] <- sd(x[i,],...)  
    }
    out <-  as.numeric(out[,1])
}
  return(out)
}






#' Variance
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{VARIANCE} Function. 
#'
#' @usage computeVariance (x,...)
#'
#' @param x Numeric or numeric vector
#' @param ... further arguments passed to or from other methods.
#'
#' @return Numeric. Returns the variance of its arguments that have valid values.
#' If the data contains missing values, it is possible to specify an na remove command, this is \code{na.rm}. The default for na.rm is \code{na.rm=F}, if the value get changed to \code{na.rm=T} every existing missing value get omitted.
#' @author Bastian Wiessner
#' @seealso \code{\link{var}}
#' @keywords internal
#' @examples
#' xpssCompute(x = fromXPSS, variables = c("V5","V7_1"),fun = "computeVariance")
#' @export

computeVariance <- function(x,...){
  if(is.matrix(x)){
    out <- x
    for(i in 1:nrow(x)){
      out[i,] <-  var(x[i,],...)  
    }
    out <-  as.numeric(out[,1])
  }
  return(out)
}
