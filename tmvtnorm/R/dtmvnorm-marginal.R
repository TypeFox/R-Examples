# Dichtefunktion und Verteilung einer multivariate truncated normal
#
# Problem ist die Bestimmung der Randverteilung einer Variablen.
#
# 1. Im bivariaten Fall kann explizit eine Formel angegeben werden (vgl. Arnold (1993))
# 2. Im multivariaten Fall kann ein Integral angegeben werden (vgl. Horrace (2005))
# 3. Bestimmung der Dichtefunktion über das Integral möglich?
# 4. Kann die Verteilungsfunktion pmvnorm() helfen? Kann man dann nach einer Variablen differenzieren?

# Literatur:
#
# Genz, A. (1992). Numerical computation of multivariate normal probabilities. Journal of Computational and Graphical Statistics, 1, 141-150
# Genz, A. (1993). Comparison of methods for the computation of multivariate normal probabilities. Computing Science and Statistics, 25, 400-405
# Horrace (2005).
# Jack Cartinhour (1990): One-dimensional marginal density functions of a truncated multivariate normal density function
# Communications in Statistics - Theory and Methods, Volume 19, Issue 1 1990 , pages 197 - 203

# Dichtefunktion für Randdichte f(xn) einer Truncated Multivariate Normal Distribution,
# vgl. Jack Cartinhour (1990) "One-dimensional marginal density functions of a truncated multivariate normal density function"
#
# @param xn Vektor der Länge l von Punkten, an dem die Randdichte ausgewertet wird
# @param i Index  (1..n) dessen Randdichte berechnet werden soll
# @param mean (nx1) Mittelwertvektor
# @param sigma (nxn)-Kovarianzmatrix
# @param lower,upper Trunkierungsvektor lower <= x <= upper
dtmvnorm.marginal <- function(xn, n=1, mean=rep(0, nrow(sigma)), sigma=diag(length(mean)),
 lower=rep(-Inf, length = length(mean)), 
 upper=rep( Inf, length = length(mean)),
 log=FALSE)
{
   if (NROW(sigma) != NCOL(sigma)) {
     stop("sigma must be a square matrix")
   }

   if (length(mean) != NROW(sigma)) {
    stop("mean and sigma have non-conforming size")
   }
   
   # Anzahl der Dimensionen                
   k <- length(mean)
   
   if (n < 1 || n > length(mean) || !is.numeric(n) || length(n) > 1 ||  !n %in% 1:length(mean))
   {
     stop("n must be a integer scalar in 1..length(mean)")
   }
   
   # Univariater Fall, vgl. Greene (2003), S.573
   if (k == 1) {
     prob    <- pnorm(upper, mean=mean, sd=sqrt(sigma)) - pnorm(lower, mean=mean, sd=sqrt(sigma))
	   density <- ifelse(
			     lower[1]<=xn & xn<=upper[1],
			     dnorm(xn, mean=mean, sd=sqrt(sigma)) / prob,
			     0)
	   if (log == TRUE) {
       return(log(density))
     } else {
	     return(density) 
	   } 
   }
   
   # Standardize sigma to correlation matrix, mean to zero vector
   # adjust xn, lower, upper
   #sd <- sqrt(diag(sigma))
   #xn    <- (xn - mean) / sd
   #lower <- (lower - mean) / sd
   #upper <- (upper - mean) / sd
   #mean  <- rep(0, k)
   #sigma  <- cov2cor(sigma)

   # Kovarianzmatrix; nach Standardisierung Korrelationsmatrix
   C <- sigma
    
   # Inverse Kovarianzmatrix, Precision matrix
   A <- solve(sigma)

   # Partitionierung von A und C
   A_1  <- A[-n,-n] # (n-1) x (n-1)
   #a_nn <- A[n, n]  # 1x1
   #a    <- A[-n, n] # (n-1) x 1
   A_1_inv <- solve(A_1)
   
   C_1  <- C[-n,-n] # (n-1) x (n-1)
   c_nn <- C[n, n]  #  1x1
   c    <- C[-n, n] # (n-1) x 1
   
   # Partitionierung von Mittelwertvektor mu
   mu   <- mean
   mu_1 <- mean[-n]
   mu_n <- mean[n]
   
   # Skalierungsfaktor der Dichte
   p <- pmvnorm(lower=lower, upper=upper, mean=mu, sigma=C)

   f_xn <- c()
   for (i in 1:length(xn))
   {
     if (!(lower[n]<=xn[i] && xn[i]<=upper[n]) || is.infinite(xn[i]))
     {
       f_xn[i] <- 0
       next
     }
     
     # m(x_n) --> (n-1x1)
     # Aufpassen bei z.B. m=c(Inf, Inf, NaN) und c=0
     m <- mu_1 + (xn[i] - mu_n) * c / c_nn
     
     # SW: Possibly optimize with vectorized version of pmvnorm() which accepts different bounds
     # for univariate density, pmvnorm() does not accept corr=
     f_xn[i] <- exp(-0.5*(xn[i]-mu_n)^2/c_nn) * pmvnorm(lower=lower[-n], upper=upper[-n], mean=m, sigma=A_1_inv)
   }
   density <- 1/p * 1/sqrt(2*pi*c_nn) * f_xn
   if (log == TRUE) {
	   return(log(density))
   } else {
	   return(density)
   }
}

