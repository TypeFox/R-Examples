jNonAdjHistVaR <-
function(s,alpha){
# Simple empirical quantile
object <- quantile(s,  alpha, names = FALSE, type = 1)
return(-object)
}
