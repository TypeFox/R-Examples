jGarchHistVaRLim <-
function(s,L, U, alpha, h){
object <- jNonAdjHistVaR(jGarchHistRetLim(s,h,L,U),alpha)
return(object)
}
