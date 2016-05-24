`alpha2rho` <-
function(alpha) {
   (alpha^2 -1 - 2* alpha * log(alpha)) / (alpha - 1) ^2
}

`rho2alpha` <-
function(rho) {
  if (abs(rho)>1) stop("rho must be within interval (0, 1)")
  f <- function(alpha, rho) {
    (rho - alpha2rho(alpha))^2
  }
  if (rho==1) {
    Inf ## must be handled in subsequent procedures
  } else if (rho==-1) {
    0
  } else if (rho==0) {
    1
  } else {
    mm <- optimize(f, interval=c(1e-12, 1e12), rho=rho)
    mm$minimum
  }
}

`pcu` <-
function(x, alpha = rho2alpha(rho), rho){
  if (is.infinite(alpha)) alpha <-  1e99  
  n <- length(x)
  y <- runif(n)
  b <- alpha+(alpha-1)*(alpha-1)*y*(1-y)
  c <- 2*y*(1-y)*(alpha*alpha*x-x+1)-2*alpha*y*(1-y)+alpha
  d <- sqrt(alpha*alpha+4*alpha*(1-alpha)*(1-alpha)*x*y*(1-x)*(1-y))
  y <- (c-(1-2*y)*d)/(2*b)
}

`pcuseries` <-
function(n, alpha = rho2alpha(rho), rho, min=0, max=1) {
   x<-numeric(n)
   x[1] <- runif(1)
   for(i in 2:n){
     x[i] <- pcu(x[i-1], alpha)
   }
   x * (max - min) + min
}

