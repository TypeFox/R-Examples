# return a vector of crps for dressed ensembles
DressCrps <- function(dressed.ens, obs) {

  if (class(obs) == "data.frame") {
    obs <- c(as.matrix(obs))
  }

  N <- nrow(dressed.ens[["ens"]])
  K.vec <- rowSums(!is.na(dressed.ens[["ens"]]))

  if (is.loaded("dresscrps")) {
    # C implementation
    crps <- with(dressed.ens, {
      sapply(1:N, function(ii) {
        .C("dresscrps", PACKAGE="SpecsVerification", 
           as.double(ens[ii,]), as.integer(K.vec[ii]),
           as.double(ker.wd[ii,]), as.double(obs[ii]), tmp=double(1)
          )[["tmp"]]
      })
    })
  } else {
    # native R implementation (slow!)
    crps <- with(dressed.ens, {
      sapply(1:N, function(ii) {
        s <- ker.wd[ii, ]
        K <- K.vec[ii]

        g1 <- data.frame(e.i=ens[ii, ], ei.y=(ens[ii, ] - obs[ii]), s.i=s)
        crps.i <- with(g1, ei.y * (2 * pnorm(ei.y / s.i) - 1)
                  + 2 * s.i * dnorm(ei.y / s.i)- s.i / K / sqrt(pi))
        g2 <- data.frame(expand.grid(i=1:K, j=1:K))
        g2 <- g2[with(g2, j < i), ]
        g2$ei.ej <- with(g2, ens[ii, j] - ens[ii, i])
        g2$si.sj <- with(g2, sqrt(s[i] ^ 2 + s[j] ^ 2))
        crps.ij <- with(g2, ei.ej * (2 * pnorm(ei.ej / si.sj) - 1)
                   + 2 * si.sj * dnorm(ei.ej / si.sj))
        sum(crps.i) / K - sum(crps.ij) / K / K
      })
    })
  }

  crps
}

