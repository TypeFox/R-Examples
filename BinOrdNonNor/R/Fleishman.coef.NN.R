Fleishman.coef.NN <-
function(skew.vec, kurto.vec)
{ 
  
  if (length(skew.vec) != length(kurto.vec)) {
    stop("Lengths of the skewness vector and the kurtosis vector differ!\n")
  }
  
  coef <- matrix(NA, length(skew.vec), 4)
  
  fleishman.poly <- function(p, r1, r2) {
    b <- p[1]
    c <- p[2]
    d <- p[3]
    r <- rep(NA, 3)
    r[1] <- b^2+6*b*d+2*c^2+15*d^2-1
    r[2] <- 2*c*(b^2+24*b*d+105*d^2+2)-r1 
    r[3] <- b*d+c^2*(1+b^2+28*b*d)+d^2*(12+48*b*d+141*c^2+225*d^2)-r2/24
    r
  }
    
  for (i in 1:length(skew.vec)) {
    r1 <- skew.vec[i]
    r2 <- kurto.vec[i]
    if (r2 < r1^2-2) stop("The skewness (r1) and kurtosis (r2) parameters for the ", i, "-th continuous variable violates the inequality: r2 >= r1^2-2!\n")
    p0 <- rep(0,3)  # starting value
    poly.coef <- BBsolve(par=p0, fn=fleishman.poly, r1=r1, r2=r2)
    if (poly.coef$conv == 0) coef[i,] <- c(-poly.coef$par[2], poly.coef$par)                          
  }
  
  # Rounding the coefficients to 9 digits after the decimal
  coef <- round(coef, 9)
  colnames(coef) <- c("a", "b", "c", "d")
  
  return(coef)
}
