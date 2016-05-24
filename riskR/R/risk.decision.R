risk.decision <-
function(x, alpha = c(0.05), beta = 1, p = 2)
{
 x = as.matrix(x)
 D = matrix(rep(0, 11 * length(alpha)), ncol = length(alpha))
 M = array(rep(0, 11 * length(alpha) * dim(x)[2]), dim = c(11, length(alpha), dim(x)[2]))
 N = rep(0, dim(x)[2])
 for (m in 1 : dim(x)[2])
 {
  M[,,m] = mean(x[,m]) / risk(x[,m], alpha, beta, p)
  N[m] = mean(x[,m])
 }
 for (j in 1 : 11)
 {
  for (k in 1 : length(alpha))
  { 
   w = which(M[j,k,] == max(M[j,k,]))
   D[j,k] = if (length(w) == 1)
   {
    w
   }else{
    which(N[w] == max(N[w]))
   } 
  }
 }
 colnames(D) <- paste(round(100*alpha, 2), "%", sep="")
 rownames(D) <- c("StD", "VaR", "EL", "ELD", "ES", "SDR", "EVaR", "DEVaR", "ENT", "DENT", "ML")
 return(D)
}
