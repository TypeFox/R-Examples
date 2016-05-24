CoxEL <- function(y, d, Z, beta, lam, fun){
#### y, d, Z  are the data inputs.(times, status and covariates)
#### beta = (beta1, ..., betak).
#### Z is matrix of covariates, size nxk; Z=(Z_1i, ..., Z_ki)
#### lam and fun controls the constrain on the baseline hazard.

yorder <- order(y, -d)
ysort <- y[yorder]
dsort <- d[yorder]

Z <- as.matrix(Z)
Zsort <- as.matrix(Z[yorder,]) ## could also use Wdataclean5( )

Zbeta <- as.vector( Zsort %*% beta )
gam <- as.vector( exp(Zbeta) )
n <- length(dsort)
fvec <- fun(ysort)

gweight <- cumsumsurv(gam)  ### rev( cumsum(rev(gam)) )  3/2015 MZ
Hw <- rep(0, n)

Hw <- 1/( gweight + lam*fvec )
Hw[ dsort == 0 ] <- 0

## for(k in 1:n) {
##  if(dsort[k] == 1) Hw[k] <- 1/( gweight[k]+lam*fvec[k] )
##  }

mu <- sum( Hw*fvec )
Hti <- cumsum(Hw)
part1 <- sum( log(Hw[Hw>0]) )
part2 <- sum( Zbeta[Hw>0] )
part3 <- - sum( gam * Hti)
part4 <- - sum( log(gweight[Hw>0]))
list(d=dsort, Hazw=Hw, mu=mu, logPlik=part2+part4, logEmpLik=part1+part2+part3)
}
