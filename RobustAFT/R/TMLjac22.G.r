TMLjac22.G <-
function(d.Beta,d.sigma,rs0,delta,X,cu) {
# Jacobian of sigma TML equation wrt sigma.hat
n     <- length(rs0); cl <- -cu
D1    <- D2 <- D3 <- D <- rep(0,n); p <- ncol(X); zero <- 1e-6
rsd   <- (rs0-X%*%d.Beta)/d.sigma
Fo    <- pnorm(rsd)
den   <- 1-Fo; ok <- den > zero
fo    <- dnorm(rsd)
ai    <- (pmax(rs0,cl)-X%*%d.Beta)/d.sigma
bi    <- (cu  -        X%*%d.Beta)/d.sigma
foai  <- dnorm(ai)
fobi  <- dnorm(bi)
Foai  <- pnorm(ai)
Fobi  <- pnorm(bi)
fopai <- foai*(-ai)
fopbi <- fobi*(-bi)
D1    <- - delta*ww(rs0,cl,cu)*2*rsd^2
D2[ok]<- -(  (1-delta)*fo/(1-Fo)^2*(foai*ai - fobi*bi + Fobi - Foai)*rsd  )[ok]
D3[ok]<- -(  (1-delta)/(1-Fo)*( fopai*ai^2 - fopbi*bi^2 )                 )[ok]
D     <- D1 + D2 + D3
Jac   <- sum(D)/d.sigma/(n-p)
Jac}

