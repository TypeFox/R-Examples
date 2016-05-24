if (!exists(".Random.seed", where = 1, mode = "numeric")) {
    rmseed <- TRUE
} else {
    saved.seed <- get(".Random.seed", pos = 1, mode = "numeric")
    rmseed <- FALSE
}
set.seed(456, kind = "Mersenne-Twister", normal.kind = "Inversion")
toy.all <- matrix(rnorm(1500*11), 1500, 11)
toy.all <- cbind(toy.all[, 2:4] %*% ((1:3)/sqrt(sum((1:3)^2))), toy.all)
rownames(toy.all) <- NULL
colnames(toy.all) <- c("y", "noise", sprintf("X%d", 1:10))
toy.learn <- toy.all[1:1000, ]
toy.test <- toy.all[1001:1500, ]
rm("toy.all")
if (rmseed) {
    rm(".Random.seed", pos = 1)
} else {
    assign(".Random.seed", saved.seed, pos = 1)
}
rm("rmseed")
