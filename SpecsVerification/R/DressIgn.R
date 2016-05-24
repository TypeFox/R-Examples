# return a vector of Ignorance scores for dressed ensembles
DressIgn <- function(dressed.ens, obs) {
  N <- nrow(dressed.ens[["ens"]])

  ign <- with(dressed.ens, {
    sapply(1:N, function(ii) {
      s <- as.numeric(ker.wd[ii, ])
      e <- as.numeric(ens[ii, ])
      o <- as.numeric(obs[ii])
      -log2(mean(dnorm(o, e, s), na.rm=TRUE))
    })
  })

  # return
  ign
}


