risk.roll <-
function(x,  N = length(x) - 1, alpha = c(0.05), beta = 1, p = 2)
{
 F = length(x) - N
 R = array(rep(0, 11 * length(alpha) * F), dim = c(11, length(alpha), F))
 for (j in 1 : F)
 {
  R[,,j] = risk(x[j : (N + j - 1)], alpha, beta, p)   
 }
 colnames(R) <- paste(round(100*alpha, 2), "%", sep="")
 rownames(R) <- c("StD", "VaR", "EL", "ELD", "ES", "SDR", "EVaR", "DEVaR", "ENT", "DENT", "ML")
 return(R)
}
