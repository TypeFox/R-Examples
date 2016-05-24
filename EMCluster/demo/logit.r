library(EMCluster, quiet = TRUE)
set.seed(1234)

x <- as.matrix(iris[, 1:4])
p <- ncol(x)
min.n <- p * (p + 1) / 2
emobj <- init.EM(x, nclass = 3, min.n = min.n, method = "Rnd.EM")

### Get posterior
post.z <- postPI(x, emobj)

### Get cov matrix of post z.
cov.param <- get.cov.param(x, emobj, post.z)
cov.post.z <- get.cov.post.z(x, emobj, post.z, cov.param = cov.param$cov)
cov.logit.z <- get.cov.logit.z(x, emobj, post.z, cov.param = cov.param$cov,
                               cov.post.z = cov.post.z)

### Get cov matrix of mixing proportion.
cov.logit.PI <- get.cov.logit.PI(x, emobj, cov.param = cov.param$cov)

### Get log or. This may have NaN because post prob can be too close to 0 or 1.
lor <- get.logor.stat(x, emobj, post.z)
