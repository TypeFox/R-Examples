jNormalVaR <-
function(s,alpha){
# Empirical estimation of normal distribution
mu <- mean(s)
vo <- sd(s)
object <- mu + vo * jNormInv(alpha)
return(-object)
}
