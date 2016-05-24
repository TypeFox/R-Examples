sizes.confint.welch <- function(alpha, ratio, delta){
  c <- ratio
  fstar <- function(nx){
    (1+c)^2/(1/(nx-1)+c^2/(c*nx-1))
  }
  g <- function(nx){
    d=delta^2
    d*nx-(1+c)*(qt(1-alpha/2,fstar(nx))^2)
  }
  if (c < 1) {a <- 1+1/c}
  if (c > 1) {a <- 2}
  if (c == 1) 
    stop("for equal variances apply 'power.t.test'")
  nx <- uniroot(g,c(a,10^6))$root
  nx <- ceiling(nx)
  ny <- c*nx
  ny <- ceiling(ny)
  structure(list(nx = nx, ny = ny, alpha = alpha, 
                 ratio = ratio, delta = delta))
}         
