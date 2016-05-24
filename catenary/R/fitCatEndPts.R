#' Get catenary parameters for given endpoints and length
#' 
#' Takes endpoints in data frame and length and returns parameters of 
#' catenary 
#'
#' @param endpoints 2x2 data frame with first column x 
#' second column is y
#' @param L length
#' @return vector of parameters
#' @author Jono Tuke, Matthew Roughan
#' @export
#' @note February 11 2013
#' @keywords internal
#' @examples
#' x <- c(-1,1)
#' y <- c(2,2)
#' endpoints <- data.frame(x=x,y=y)
#' L <- 3
#' fitCatEndPts(endpoints,L)
fitCatEndPts <- function(endpoints,L){
  if(L < minmaxLength(endpoints)){
    stop("Not possible: length is too short")
  }
  fn <- function(para,endpoints,L){
    x0 <- endpoints[1,1];x1 <- endpoints[2,1]
    y0 <- endpoints[1,2];y1 <- endpoints[2,2]
    c1 <- para[1];c2 <- para[2]; lambda <- para[3]
    
    d1 <- (y0 - f(x=x0,c1=c1,c2=c2,lambda=lambda))^2
    d2 <- (y1 - f(x=x1,c1=c1,c2=c2,lambda=lambda))^2
    d3 <- (L - getCatLength(x0=x0,x1=x1,c1=c1,c2=c2))^2
    
    return( sum(c(d1,d2,d3)))
  }
  c2 <- mean(endpoints[,1])
  c1 <- min(endpoints[,2])/2
  lambda <- c1
  
  tryCatch(tmp <- optim(par=c(c1,c2,lambda),fn=fn,endpoints=endpoints,L=L),
           error = function(e) e,
           finally=print("optim worked"))
  if(tmp$convergence !=0){
    print("Maybe problems with convergence of optim - be careful")
  }
  par <- c(c1=tmp$par[1],c2=tmp$par[2],lambda=tmp$par[3])
  return(par)
}
