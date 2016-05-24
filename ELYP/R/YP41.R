YP41 <- function(y, d, Z, b1, b2, k){
#### y, d, Z  are the data inputs.(times, status and covariates)
#### b1, b2 needs to be given.
#### b1 = (beta1) and b2=(beta2).
#### Z is vector of covariates, size nx1; Z=(Z_i)
#### k(+1) is the number of iterations we do, in maximizing the baseline.
#### for model without alpha

Z <- as.vector(Z)
yorder <- order(y, -d)
ysort <- y[yorder]
dsort <- d[yorder]
Zsort <- Z[yorder]  ## could also use Wdataclean5( )

n <- length(Zsort)        #### length of y, d should be also n
S0 <- (n:1)/(n+1)   #### is there a better initial value?

gam <- exp( - as.matrix(Zsort) %*% cbind(b1,b2) )
tempT <- Haz4(d=dsort, S=S0, gam=gam)

#tempT <- Haz4(d=dsort, S=temp3$Su, gam=gam)

i <- 1
err <- 8
while( i <= k & err > 0.0000001 ) {
OldSu <- tempT$Su
gam <- exp( - as.matrix(Zsort) %*% cbind(b1, b2) )
tempT <- Haz4(d=dsort, S=tempT$Su, gam=gam)
i <- i+1
err <- sum( abs(tempT$Su - OldSu) )
}

if(i == (k+1)) warning("check convergence for Baseline S(t)")

list(d=dsort, Hazw=tempT$Hazw, Survival=tempT$Su, gam=gam)
}
