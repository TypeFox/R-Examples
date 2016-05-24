Y <- rnorm(20L)
Z <- rnorm(20L)
cor(Y, Z, method = "spearman")

ranksY <- rank(Y)
ranksZ <- rank(Z)
cor(ranksY, ranksZ, method = "pearson")

