#' Get catenary parameters for given endpoints and length
#' 
#' Takes endpoints in data frame and length and returns parameters of 
#' catenary 
#'
#' @param endpoints 2x2 data frame with first column x 
#' second column is y
#' @return vector of parameters
#' @author Jono Tuke <simon.tuke@@adelaide.edu.au>
#' @export
#' @note February 11 2013
#' @keywords internal
#' @examples
#' x <- c(-1,1)
#' y <- c(2,2)
#' endpoints <- data.frame(x=x,y=y)
#' fitNaturalCat(endpoints)
fitNaturalCat <- function(endpoints){
  c1 <- min(endpoints[,2])
  c2 <- mean(endpoints[,1])
  
  fn <- function(para,endpoints){
    x0 <- endpoints[1,1];x1 <- endpoints[2,1]
    y0 <- endpoints[1,2];y1 <- endpoints[2,2]
    c1 <- para[1];c2 <- para[2]; lambda <- 0
    
    d1 <- (y0 - f(x=x0,c1=c1,c2=c2,lambda=lambda))^2
    d2 <- (y1 - f(x=x1,c1=c1,c2=c2,lambda=lambda))^2
    
    return( sum(c(d1,d2)))
  }
  tryCatch({
    tmp <- optim(par=c(c1,c2),fn=fn,endpoints=endpoints)
  }, error = function(e){
    print(e)
  }, finally= {
    cat("Optim worked\n")
  })
  if(tmp$convergence !=0){
    print("Maybe problems with convergence of optim - be careful")
  }
  par <- c(c1=tmp$par[1],c2=tmp$par[2],lambda=0)
  return(par)
}
