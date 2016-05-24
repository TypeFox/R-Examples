jGarch11VaRLim <-
function(s,L,U,alpha,h){
# Garch(1,1)
x0 <- c(0.05,0.85,mean(s^2)*0.1)
 x <- jNewRapGarchLim(x0, s, L, U)
object <- sqrt(jGarch11VolLim(s,x,h)) * jNormInv(alpha)
return(-object)
}
