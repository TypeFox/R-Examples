mplik.wb.bi <-
function(par,Y,X,delta,whc) {
          y <- log(Y)
      sigma <- exp(par[1])
        phi <- matrix(par[-1],ncol=1)
          r <- sum(delta)
          z <- (y-X%*%phi)/sigma
loglikelihood <- -r*log(sigma)+sum((y-X%*%phi)*delta)/sigma-sum(exp((y-X%*%phi)/sigma))
J <- J.inf.weibul(Y, X, sigma, phi, delta, whc)
LX <- LX.mat.weibull(Y, X, sigma, phi, delta, whc)
mplk <- loglikelihood + .5*log(abs(det(J))) -log(abs(det(LX)))
return(-mplk)
}
