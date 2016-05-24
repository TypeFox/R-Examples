# Verteilungsfunktion für die eindimensionale Randdichte f(xn) einer Truncated Multivariate Normal Distribution,
# vgl. Jack Cartinhour (1990) "One-dimensional marginal density functions of a truncated multivariate normal density function" für die Dichtefunktion
#
# @param xn Vektor der Länge l von Punkten, an dem die Verteilungsfunktion ausgewertet wird
# @param i Index  (1..n) dessen Randdichte berechnet werden soll
# @param mean (nx1) Mittelwertvektor
# @param sigma (nxn)-Kovarianzmatrix
# @param lower,upper Trunkierungsvektor lower <= x <= upper
ptmvnorm.marginal <- function(xn, n=1, mean=rep(0, nrow(sigma)), sigma=diag(length(mean)), lower=rep(-Inf, length = length(mean)), upper=rep( Inf, length = length(mean)))
{
   # check of standard tmvnorm arguments
   cargs <- checkTmvArgs(mean, sigma, lower, upper)
   mean  <- cargs$mean
   sigma <- cargs$sigma
   lower <- cargs$lower
   upper <- cargs$upper
   
   if (n < 1 || n > length(mean) || !is.numeric(n) || length(n) > 1 ||  !n %in% 1:length(mean))
   {
     stop("n must be a integer scalar in 1..length(mean)")
   }
   
   # Anzahl der Dimensionen                
   k = length(mean)
   
   Fx     = numeric(length(xn))
   upper2 = upper
   alpha  = pmvnorm(lower = lower, upper = upper, mean = mean, sigma = sigma)
   for (i in 1:length(xn))
   {
    upper2[n] = xn[i]
    Fx[i]     = pmvnorm(lower=lower, upper=upper2, mean=mean, sigma=sigma)
   }
   return (Fx/alpha)
}
