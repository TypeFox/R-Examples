#' Bidimentional Sinc train set example
#'
#' Generates the training set of sinc(x)*sinc(y) for the (x,y) regular grid
#'
#' 
#' @param x numeric vector with the x-th grid coordinates
#' @param y numeric vector with the x-th grid coordinates
#'
#' 
#' @return
#' \item{matrix}{numeric matrix with the columns x, y and z=sync(x,y)} 
#'
#' @include Anfis-training.R
#' @export
#' @name trainSet
#' @rdname ANFIS-trainSet
#' @family ANFIS
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
#' @examples
#' ##Domain definition for a regular (x,y) grid with 11 points for each 
#' ##coordinates
#' x <- seq(-10, 10, length= 11)
#' trainingSet <- trainSet(x,x)
#' Z <- matrix(trainingSet[,"z"],ncol=length(x),nrow=length(x))
#' 
#' ##Plot the domain
#' persp(x, x, Z, theta=45, phi=15, expand=0.8, col="lightblue", 
#'  ticktype="detailed", main="sinc(x)*sinc(y)")
trainSet <- function(x, y){
  sinc <- function(z){
    ifelse(z==0,1,sin(z)/z)
  }
  out <- matrix(ncol=3,nrow=0)
  colnames(out) <- c("x", "y", "z")
  invisible(sapply(x, function(X){
    sapply(y, function(Y){
      out <<- rbind(out, c(X, Y, sinc(X)*sinc(Y)))
      return(NULL)
    })
  }))
  return(out)
}
