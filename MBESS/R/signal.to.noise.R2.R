signal.to.noise.R2 <- function(R.Square, N, p)
{
Y <- Y.Adj <- Theta.U <- Theta.L <- Theta.NL <- NA

# From Muirhead, 1985
# Uses his notation but define his parameters in terms of N and p.
n <- N-1
m <- p+1

Y <- R.Square/(1-R.Square)

# Often used method of adjustment.
R.Square.Adj <- ((N-1)/(N-p-1))*R.Square - (p)/(N-p-1)

Y.Adj <- R.Square.Adj/(1-R.Square.Adj)

if(n-m-3 >=1)
{
# Theta.U is equal to Stuart, Ord & Arnold's (1999) Equation 28.97
Theta.U <- ((n-m-1)/n)*Y - (m-1)/n
Theta.U <- max(Theta.U, 0)

# Theta.L is equal to Stuart, Ord & Arnold's (1999) Equation 28.98
Theta.L <- ((n*(n-m-3))/((n+2)*(n-m-1)))*Theta.U
Theta.L <- max(Theta.L, 0)

# Only for m >= 6.
if(m>=6)
{
Theta.NL <- Theta.L + ((2*(n-2)*(n-m-3)*(m-5))/((n+2)*(n-m-1)*(n-m+1)*(n-m+3)))*(1/Y)
Theta.NL <- max(Theta.NL, 0)
}
else(Theta.NL <- "p must be >= 5 in order to calculate \'phi2.UMVUE.NL\'")
}
return(list(phi2.hat=Y, phi2.adj.hat=Y.Adj, phi2.UMVUE=Theta.U, phi2.UMVUE.L=Theta.L, phi2.UMVUE.NL=Theta.NL))
}
