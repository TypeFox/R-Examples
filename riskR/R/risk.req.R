risk.req <-
function(x, M = 10^6, T = 1, alpha = c(0.05), beta = 1, p = 2)
{
 R <- risk(x, alpha, beta, p) * M * sqrt(T)
 colnames(R) <- paste(round(100*alpha, 2), "%", sep="")
 rownames(R) <- c("StD", "VaR", "EL", "ELD", "ES", "SDR", "EVaR", "DEVaR", "ENT", "DENT", "ML")
 return(R)
}
