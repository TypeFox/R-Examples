jNormInv <- 
function(u){
# Absolute error .45e-3  , Hastings (1955).
c0 <- 2.515517
c1 <- .802853
c2 <- .010328
d1 <- 1.432788
d2 <- .189269
d3 <- .001308
z <- (-2*log(pmin(u,1-u)))^.5
object <- sign(u-0.5)*(z - (c0 + c1*z + c2*z^2)/(1 + d1*z + d2*z^2 + d3*z^3))
return(object)
}
