#' @title Multi-collinearity test
#' @description Test for multi-collinearity in data using qr-matrix decomposition
#'
#' @param x data.frame or matrix object 
#' @param p multi-collinearity threshold
#'
#' @return test statistic message
#' @return names character vector of multi-collinearity variables
#'
#' @author Jeffrey S. Evans  <jeffrey_evans<at>tnc.org>
#'
#' @references
#'  Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S Language. Wadsworth & Brooks/Cole. 
#' @references
#'  Dongarra, J. J., Bunch, J. R., Moler, C. B. and Stewart, G. W. (1978) LINPACK Users Guide. Philadelphia: SIAM Publications. 
#'
#' @note
#' The multi-collinearity threshold needs to be adjusted based on number of parameters. For small number(s) of variables (<20) use ~1e-07 and for larger ~0.05  
#'
#' @examples  
#' test <- data.frame(v1=seq(0.1, 5, length=100), v2=seq(0.1, 5, length=100), 
#'                    v3=dnorm(runif(100)), v4=dnorm(runif(100))) 
#'   ( cl <- multi.collinear(test) )
#'
#' # PCA biplot of variables
#' pca.test <- prcomp(test[,1:ncol(test)], scale=TRUE)
#'   biplot(pca.test, arrow.len=0.1, xlabs=rep(".", length(pca.test$x[,1])))        
#'
#'  # Remove identified variable(s)
#'  test <- test[,-which(names(test)==cl)]
#'
#' @export
multi.collinear <- function(x, p=1e-07) {
  if (!inherits(x, "data.frame")) stop("X MUST BE A data.frame")
  if ( (dim(x)[2] < 2) == TRUE) stop("NOT ENOUGH VARIABLES TO TEST")
    xtest <- x
    x.names <- names(xtest)
    qrx <- qr(xtest, tol=p)
    if (length(names(xtest)[qrx$pivot[1:qrx$rank]]) != length(xtest) )  
      {  
        keep <- names(xtest)[qrx$pivot[1:qrx$rank]]
         warning("MULTI-COLINEAR VARIABLES: ", paste(setdiff(x.names, keep), collapse = ","))
      return(paste(setdiff(x.names, keep)))
    } else { print(" NO MULTICOLINEAR VARIABLES IDENTIFIED")
  } 
}
