#' @title Compare a sample to one or more fitted distributions
#'
#' @param X An unweighted sample
#' @param Dist1,Dist2,Dist3 The fitted distribution, specified as either the  objects of class eDist or names of the distribution to
#' be fitted.
#' @author Haizhen Wu and A. Jonathan R. Godfrey.
#'  
#' @return compareDist returns an object of class histogram comparing the sample distribution to the specified fitted distribution(s).
#' @examples
#' X <- rBeta(n=100, params=list(shape1=1, shape2=2))
#' compareDist(X, "Beta", "Normal", eNormal(X))


#' @export compareDist
compareDist <- function(X, Dist1, Dist2=NULL, Dist3=NULL){
  AllowDists <- c("Normal","Beta") # add in the other distributions here, perhaps using an internal object
  if(class(Dist1)!="eDist"){
    if(is.element(Dist1, AllowDists)){ 
      # get the eNormal etc details here
      eFoo <- get(paste("e", Dist1, sep=""))
      Dist1 <- eFoo(X)
    }
    else stop("No matching distribution name found for Dist1.\n")
  }
  
  if(!is.null(Dist2)){
    if(class(Dist2)!="eDist"){
      if(is.element(Dist2, AllowDists)){ # add in the other distributions here
        # get the eNormal etc details here
        eFoo <- get(paste("e", Dist2, sep=""))
        Dist2 <- eFoo(X)
      }
      else stop("No matching distribution name found for Dist2.\n")
    }
  }
  if(!is.null(Dist3)){
    if(class(Dist3)!="eDist"){
      if(is.element(Dist3, AllowDists)){ # add in the other distributions here
        # get the eNormal etc details here
        eFoo <- get(paste("e", Dist3, sep=""))
        Dist3 <- eFoo(X)
      }
      else stop("No matching distribution name found for Dist3.\n")
    }
  }
  
  n <- length(X)
  xmin <- min(X)
  xmax <- max(X)
  
  if(n>50){
    # do histogram
    hist(X, freq=FALSE, main="")
    dFoo <- function(x) {get(paste("d", attributes(Dist1)$distname, sep=""))(x, params=Dist1)}
    curve(dFoo, from=xmin, to=xmax, add=T, col=1)
    if(!is.null(Dist2)){
      dFoo <- function(x) {get(paste("d", attributes(Dist2)$distname, sep=""))(x, params=Dist2)}
      curve(dFoo, from=xmin, to=xmax, add=T, col=2, lty=2)
      if(!is.null(Dist3)){
        dFoo <- function(x) {get(paste("d", attributes(Dist3)$distname, sep=""))(x, params=Dist3)}
        curve(dFoo, from=xmin, to=xmax, add=T, col=4, lty=3)
      }
    }
    
  }
  else{
    # do q-q or p-p plot
  }
  
}
