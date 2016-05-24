#' Get catenary parameters for given endpoints and length
#' 
#' Takes endpoints in data frame and length and returns parameters of 
#' catenary 
#'
#' @param endpoints 2x2 data frame with first column x 
#' second column is y
#' @param L length
#' @return vector of parameters
#' @keywords internal
#' @author Jono Tuke, Matthew Roughan
#' @export
#' @note February 11 2013
#' @examples
#' x <- c(-1,1)
#' y <- c(2,2)
#' endpoints <- data.frame(x=x,y=y)
#' L <- 3
#' fitCat(endpoints,L)
fitCat <- function(endpoints,L){
  if(L < minmaxLength(endpoints)[1]){
    stop("Not possible: length is too short")
  }
  x0 <- endpoints[1,1];x1 <- endpoints[2,1]
  y0 <- endpoints[1,2];y1 <- endpoints[2,2]
  k <- sqrt(L^2 - (y1-y0)^2) / (x1-x0)
  
  fn <- function(phi,k){
    return(sinh(phi)-k*phi)
  }
  
  tryCatch({
    tmp <- uniroot(fn,interval=c(acosh(k),sqrt(6*(k-1))),k=k)
  }, error = function(e){
    stop(e)
  }, finally = {
    cat("uniroot worked\n")
  }
           )
  phi <- tmp$root
  c1 <- (x1 - x0) / (2*phi)
  theta <- atanh( (y1 - y0) / L )
  c2 <- 0.5*(x1 + x0 - theta*(x1 - x0)/phi)
  lambda1 <- y1  - c1*cosh((x1-c2)/c1)
  lambda0 <- y0  - c1*cosh((x0-c2)/c1)
  if(lambda0 != lambda1){
    stop("lambda0 and lambda1 disagrees so summ'it is wrong")
  }
  return(c(c1=c1,c2=c2,lambda=lambda0))
}
