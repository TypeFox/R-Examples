TMLjac21.W <-
function(d.Beta,d.sigma,rs0,delta,X,cl,cu) {
# Jacobian of sigma-equation wrt Beta.hat
n     <- length(rs0); zero <- 1e-6
D1    <- D2 <- D3 <- D <- rep(0,n); p <- ncol(X)
rsd   <- (rs0-X%*%d.Beta)/d.sigma
Fo    <- plweibul(rsd)
fo    <- dlweibul(rsd)
ok    <- (1-Fo) > zero
ai    <- (pmax(rs0,cl)-X%*%d.Beta)/d.sigma
bi    <- (cu  -        X%*%d.Beta)/d.sigma
foai  <- dlweibul(ai)
fobi  <- dlweibul(bi)
Foai  <- plweibul(ai)
Fobi  <- plweibul(bi)
fopai <- foai*(-ps0W(ai))
fopbi <- fobi*(-ps0W(bi))
D1    <- - delta*ww(rs0,cl,cu)*psp1W(rsd)
D2[ok]    <- - (  (1-delta)*fo/(1-Fo)^2*(foai*ai-fobi*bi + Fobi - Foai)  )[ok]
D3[ok]    <- - (  (1-delta)/(1-Fo)*( fopai*ai - fopbi*bi )               )[ok]
D     <- D1 + D2 + D3
Jac   <- t(X)%*%(as.vector(D))/d.sigma/(n-p)
Jac}

