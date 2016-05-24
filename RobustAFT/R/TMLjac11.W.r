TMLjac11.W <-
function(d.Beta,d.sigma,rs0,delta,X,cl,cu) {
# Jacobian of Beta-equation wrt Beta.hat
n     <- length(rs0)
D1    <- D2 <- D3 <- D <- rep(0,n); zero <- 1e-6
rsd   <- (rs0-X%*%d.Beta)/d.sigma
Fo    <- plweibul(rsd)
fo    <- dlweibul(rsd)
ok    <- (1-Fo) > zero
ai    <- (pmax(rs0,cl)-X%*%d.Beta)/d.sigma
bi    <- (cu  -        X%*%d.Beta)/d.sigma
foai  <- dlweibul(ai)
fobi  <- dlweibul(bi)
fopai <- foai*(-ps0W(ai))
fopbi <- fobi*(-ps0W(bi))
D1    <- - delta*ww(rs0,cl,cu)*psp0W(rsd)
D2[ok]    <- - (  (1-delta)*fo/(1-Fo)^2*(foai-fobi )  )[ok]
D3[ok]    <- - (  (1-delta)/(1-Fo)*( fopai-fopbi )    )[ok]
D     <- D1 + D2 + D3
Jac   <- t(X) %*% (as.vector(D)*X) /d.sigma/n
Jac}

