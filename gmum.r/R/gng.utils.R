#' @importFrom httr GET content

#' @export
#' @rdname print.gng
#' @method print Rcpp_GNGServer
#'
#' @title print
#' 
#' @description Print basic information about GNG object
#'
#' @docType methods
#'
#' @param x GNG object model.
#' @param ... other arguments not used by this method.
print.Rcpp_GNGServer <- NULL

#' Summary of GNG object
#' @export
#' @rdname summary.gng
#' @method summary Rcpp_GNGServer
#'
#' @title summary
#' 
#' @description Print basic information about GNG object
#'
#' @docType methods
#'
#' @param object GNG object model.
#' @param ... other arguments not used by this method.
summary.Rcpp_GNGServer <- NULL


print.Rcpp_GNGServer <- function(x, ...){
  print(sprintf("Growing Neural Gas, %d nodes with mean error %f", 
                x$getNumberNodes(), x$getMeanError()))
}

summary.Rcpp_GNGServer <- function(object, ...){
  if(object$.getConfiguration()$.uniformgrid_optimization){
    print("(Optimized) Growing Neural Gas")
  }else{
    print("Growing Neural Gas")
  }
  if(exists("object$call")){
    print(object$call)
  }
  if(object$hasStarted()){
    print(sprintf("%d nodes with mean error %f", 
                  object$getNumberNodes(), object$getMeanError()))
    
    print(sprintf("Trained %d iterations", object$getCurrentIteration()))
    print("Mean errors[s]: ")
    errors = object$getErrorStatistics()
    if(length(errors) > 10){
      errors = errors[(length(errors)-10):length(errors)]
    }
    
    print(errors)
  }
}

show.Rcpp_GNGServer <- function(object) {
  summary(object)
}

setMethod("show", "Rcpp_GNGServer", show.Rcpp_GNGServer)

#' Retrieves wine dataset design matrix from UCI repository
#' 
#' @title get.wine.dataset.X
#' 
#' @param scale if TRUE will perform feature scaling
#' 
#' @export
get.wine.dataset.X <- function(scale=TRUE){
  if(!exists(".wine") || is.null(.wine)) {
    a <- GET("https://archive.ics.uci.edu/ml/machine-learning-databases/wine/wine.data")
    .wine <<- read.csv(textConnection(content(a)), header=F)
  }
  
  if(scale) {
    return(as.matrix(scale(.wine[-1])))
  } else {
    return(.wine[-1])
  }
}

#' Retrieves wine dataset labels from UCI repository
#' 
#' @title get.wine.dataset.y
#' 
#' @export
get.wine.dataset.y <- function(){
  # Hack for R CMD check. Note that it is cleaner to assign (see predictComponent)
  if(!exists(".wine") || is.null(.wine)) {
    a <- GET("https://archive.ics.uci.edu/ml/machine-learning-databases/wine/wine.data")
    .wine <<- read.csv(textConnection(content(a)), header=F)
  }
  return(.wine[,1])
}

.plane.point<-function(r,center){
  if(!hasArg(r)) r<-1.0
  if(!hasArg(center)) center<-c(0,0,0)
  
  point<-center
  point[1]<-point[1]+r*runif(1.0)
  point[2]<-point[2]+r*runif(1.0)
  point[3]<-point[3]
  
  return(point)
}

.sphere.point<-function(r,center){
  if(!hasArg(r)) r<-1.0
  if(!hasArg(center)) center<-c(0,0,0)
  
  alpha<-runif(1)*2*pi
  beta<-runif(1)*pi
  
  point<-center
  point[1]<-point[1]+r*cos(alpha)*sin(beta)
  point[2]<-point[2]+r*sin(alpha)*sin(beta)
  point[3]<-point[3]+r*cos(beta)
  
  return(point)
}
