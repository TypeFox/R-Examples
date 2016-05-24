risk.hedge <-
function(x, y, alpha = c(0.05), beta = 1, p = 2)
{
 H = matrix(rep(0, 11 * length(alpha)), ncol = length(alpha))
 for (j in 1 : 11)
 {
  for (k in 1 : length(alpha))
  { 
   f <- function (h)
   {
    risk(x - h * y, alpha, beta, p)[j, k]
   }
   H[j,k] = optimize(f, c(-20, 20)) $ minimum
  }
 }
 colnames(H) <- paste(round(100*alpha, 2), "%", sep="")
 rownames(H) <- c("StD", "VaR", "EL", "ELD", "ES", "SDR", "EVaR", "DEVaR", "ENT", "DENT", "ML")
 return(H)
}
