jMleNormalVaR <-
function(s,alpha){
mu <- mean(s)
vo <- sd(s) * sqrt((length(s)-1)/length(s))
object <- mu + vo * jNormInv(alpha)
return(-object)
}
