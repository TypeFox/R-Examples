FtoTheta <-
function(p,f) {
   num <- 4*p*(1-p)*(1-f)^2
   den <- (p+f*(1-p))*((1-p)+f*p)
   theta <- num/den
   return(theta)
}
