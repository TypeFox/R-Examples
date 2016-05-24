gedlogl <-
function(z, p = 2)
{
lz <- length(z)
abszp <- abs(z)^p
G3p <- gamma(I(3/p))
logG3p <- lgamma(I(3/p))
G1p <- gamma(I(1/p))
logG1p <- lgamma(I(1/p))

#term 1:
term1 <- lz*(log(p) + 0.5*logG3p - log(2) - (3/2)*logG1p)

#term 2:
term2 <- sum(abszp)*(G3p/G1p)^(p/2)

#return logl:
return(as.vector(term1 - term2))
}
