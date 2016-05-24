
 library(mcmc)

 set.seed(42)

 n <- 1e5
 rho <- 0.99

 x <- arima.sim(model = list(ar = rho), n = n)
 gamma <- acf(x, lag.max = 1999, type = "covariance",
     plot = FALSE)$acf
 k <- seq(along = gamma)
 Gamma <- gamma[k %% 2 == 1] + gamma[k %% 2 == 0]
 k <- min(seq(along = Gamma)[Gamma < 0])
 Gamma <- Gamma[1:k]
 Gamma[k] < 0
 Gamma[k] <- 0

 out <- .Call("initseq", x - mean(x))
 names(out)

 all.equal(gamma[1], out$gamma0)

 length(out$Gamma.pos) == length(Gamma)
 all.equal(out$Gamma.pos, Gamma)

 Gamma.dec <- cummin(Gamma)
 all.equal(out$Gamma.dec, Gamma.dec)
 
library(Iso)
 Gamma.con <- Gamma.dec[1] + cumsum(c(0, pava(diff(Gamma.dec))))
 all.equal(out$Gamma.con, Gamma.con)

 all.equal(0, min(out$Gamma.pos - out$Gamma.dec))
 max(diff(out$Gamma.dec)) < sqrt(.Machine$double.eps)

 all.equal(0, min(out$Gamma.dec - out$Gamma.con))
 min(diff(diff(out$Gamma.con))) > (- sqrt(.Machine$double.eps))

 all.equal(2 * sum(out$Gamma.pos) - out$gamma0, out$var.pos)
 all.equal(2 * sum(out$Gamma.dec) - out$gamma0, out$var.dec)
 all.equal(2 * sum(out$Gamma.con) - out$gamma0, out$var.con)

 rev(out$Gamma.pos)[1] == 0
 rev(out$Gamma.dec)[1] == 0
 all.equal(rev(out$Gamma.con)[1], 0)

