size.variance.confint <- function(alpha, delta = NULL, deltarel = NULL){
  c1 <- function(n){
    qchisq(1-alpha/2,n-1)
  }
  c2 <- function(n){
    qchisq(alpha/2,n-1)
  }
  if (sum(sapply(list(delta,deltarel), is.null)) != 1) 
    stop("exactly one of 'delta' and 'deltarel' must be NULL")
  if(is.null(deltarel)){
    g1 <- function(n)
       {
         (n-1)*(c1(n)-c2(n))-2*delta*c1(n)*c2(n)
       }
     n <- uniroot(g1,c(2,10^6))$root
     n <- ceiling(n)
     n
   }
  if(is.null(delta)){
    g2 <- function(n)
      {
        deltarel*(c1(n)+c2(n))-(c1(n)-c2(n))
      }
    n <- uniroot(g2,c(2,10^6))$root
    n <- ceiling(n)
    n
   }
  structure(list(n = n, alpha = alpha, delta = delta, deltarel = deltarel))
}         
