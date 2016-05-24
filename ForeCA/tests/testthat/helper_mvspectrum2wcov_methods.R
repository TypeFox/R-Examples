
kNumVariables <- 3
kNumObs <- 1e4 - 1
kSeries <- matrix(rnorm(kNumObs * kNumVariables), ncol = kNumVariables)
kSeriesCentered <- sweep(kSeries, 2, colMeans(kSeries), "-")
kWhitenedSeries <- whiten(kSeries)$U

Sigma <- cov(kSeries)

ww.tmp <- initialize_weightvector(num.series = ncol(kSeries),
                                  method = "rnorm")
yy <- c(kWhitenedSeries %*% t(ww.tmp))
ec.tmp <- list(prior.weight = 0.1)
