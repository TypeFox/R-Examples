jMonteCarloVarLim <-
function(s,L,U, alpha){
# T is number of day of sim
mu <- mean(s)
vo <- sd(s)
m <- 1000
ret <- jMCPriLim(1,mu,vo,m)-1
object <- jNonAdjHistVaR(ret,alpha)
return(object)
}
