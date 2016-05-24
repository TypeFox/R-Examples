jGarch11VaR <-
function(s,alpha,h){
# Garch(1,1)
x0 <- c(0.05,0.85,mean(s^2)*0.1)
 x <- jNewRapGarch(x0, s)
object <- sqrt(jGarch11Vol(s,x,h)) * jNormInv(alpha)
return(-object)
}
