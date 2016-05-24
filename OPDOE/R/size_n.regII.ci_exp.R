size_n.regII.ci_exp <- function(mx, my, sx2=((xmax-xmin)/6)^2, sy2, xmin, xmax, alpha=0.05){
  if(xmin>xmax)
    stop("xmax < xmin !")
  if(sx2<0)
    stop("sx2<0 !")
  if(sy2<0)
    stop("sy2<0 !")

  n.3 <- 1003
  n.2 <- 1002 
  n.1 <- Inf
  # until n stays constant or jitters between two values:
  while(n.2!=n.1 & n.1!=n.3){
    n <- ceiling(qt(1-alpha/2, n.1-2)^2*(1+max(c((xmin-mx)^2, (xmax-mx)^2))/sx2))
    n.3 <- n.2
    n.2 <- n.1
    n.1 <- n
  }
  n
}
