if (!exists(".Random.seed", where = 1, mode = "numeric")) {
    rmseed <- TRUE
} else {
    saved.seed <- get(".Random.seed", pos = 1, mode = "numeric")
    rmseed <- FALSE
}
set.seed(123, kind = "Mersenne-Twister", normal.kind = "Inversion")
## Special initialization for 1st and 2nd elements:
## match autocovariance of the rest of the series
toyTs <- c(drop(matrix(c(1, -5/7, 0, sqrt(24)/7), 2, 2) %*% rnorm(2)),
           rnorm(3000, sd = sqrt(0.66 - 3 / 14)))
## Autocovariances 1, -5/7 (-0.71), 23/35 (0.66), -19/35 (-0.54),
## 82/175 (0.47), -139/350 (-0.40), 1187/3500 (0.34), -2021/7000
## (-0.29), 2461/10000 (0.25), ...
toyTs.coeffs <- c(-0.5, 0.3)
for (i in 3:3002) {
    toyTs[i] <- toyTs[i] +
        toyTs.coeffs[1] * toyTs[i-1] + toyTs.coeffs[2] * toyTs[i-2]
}
## Drop 1st and 2nd element, separate learning and test sets
tsToy.learn <- toyTs[3:1002]
tsToy.test <- toyTs[1003:3002]
rm("toyTs.coeffs", "toyTs", "i")
if (rmseed) {
    rm(".Random.seed", pos = 1)
} else {
    assign(".Random.seed", saved.seed, pos = 1)
}
rm("rmseed")
