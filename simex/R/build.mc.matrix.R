build.mc.matrix <-
function (mc.matrix, method = "series" , tuning = sqrt(.Machine$double.eps),
diag.cor = FALSE, tol = .Machine$double.eps, max.iter = 100)
{
# mc.matrix is a matrix which should exist for mc.matrix^lambda but does not
# tuning should be a very small amaount which is added due to numerical
# reasons in the log criterion
# diag.cor schould the amount of corrected values be substracted form the
# diagonal value (= TRUE) or from all positive values proportional to their
# size (= FALSE)
# tol: gives the minimum for the progression described in equation(2.1) for
# convergence
# checking if the matrix is already a valid matrix
# if (check.mc.matrix(list(mc.matrix))) stop("Matrix is already a valid mc.matrix")
nam <- dimnames(mc.matrix)
if (method == "jlt") {
if (any(diag(mc.matrix == 0)))
stop("0 or 1 on the diagonal are not allowed for method == jlt")
   # see israel et al. page 251 for details
   generator <- t(t(mc.matrix) * log(diag(mc.matrix)) / (diag(mc.matrix) - 1))
   diag(generator) <- log(diag(mc.matrix))
}
# The serie (P-I) - (P^2-I)/2 + (P^3-I)/3 + ... = Q should yield a valid
# generator for Pi and therefore a valid misclassification matrix
if (method == "series") {
iter <- 1
eps <-1
ident <- diag(dim(mc.matrix)[1]) # Identity matrix
generator <- matrix(0, nrow = dim(mc.matrix)[1], ncol = dim(mc.matrix)[2])
while (sum(abs(eps)) > tol && iter < max.iter) {
eps <- ident
for (i in 1:iter) eps <- eps %*% (mc.matrix - ident)
eps <- (-1)^(iter + 1) * eps / iter
generator <- generator + eps
iter <- iter + 1
}
if (iter == max.iter)
stop("Series did not converge, try method == jlt instead")
}
# The "log" method builds an generator voa the log which is not always valid
# but can be made valid via adding all negativ non diagonal values up to 0
if (method == "log") {
# building the log(mc.matrix)
ev <- eigen(mc.matrix)
evalue <- ev[["values"]]
evectors <- ev[["vectors"]]
d <- diag(log(evalue))
generator <- evectors %*% d %*% solve(evectors)
}
# checking each colum of log(mc.matrix) if there are negative values on
# the off-diagonal. Correcting negative values on the off-doagonal and
# setting them to 0 + tuning
# correcting the diagonal value and setting it to old value + the sums of
# the corrected values for this column (equation 3.1)
if (diag.cor == TRUE) {
for (i in 1:dim(mc.matrix)[1]) {
if (any(generator[-i, i] < 0)) {
index <-  generator[, i] < 0
corrector <- rep(0, dim(generator)[1])
corrector[index] <- abs(generator[index, i]) + tuning
index[i] <- F
corrector[i] <- 0
generator[, i] <- generator[, i] + corrector
generator[i, i] <- generator[i, i] - sum(corrector)
}
}
} else {
# correcting as above but distributing the amount of correction
# on all values proportional to their size (equation 3.1´)
index2 <- seq(1, dim(mc.matrix)[2])
for (i in 1:dim(mc.matrix)[1]) {
corrector <- rep(0, dim(generator)[1])
g <- abs(generator[i, i]) + sum(apply(cbind(generator[-i, i], 0), 1, max))
b <- sum(apply(cbind(-generator[-i, i] + tuning, 0), 1, max))
index <-  generator[, i] < 0
corrector[index] <- tuning
corrector[i] <- generator[i, i]
index2 <- g > 0
index[i] <- FALSE
corrector[!index * index2] <- (generator[, i] -
(b * abs(generator[, i]) / g))[!index * index2]
if (g == 0) corrector <- generator
generator[, i] <- corrector
}
}
# building the corrected mc.matrix
ev <- eigen(generator)
evalue <- ev[["values"]]
evectors <- ev[["vectors"]]
d <- diag(exp(evalue))
mc.matrix <- zapsmall(evectors %*% d %*% solve(evectors))
# mc.matrix[mc.matrix < 0 ] <- 0
dimnames(mc.matrix ) <- nam
return(mc.matrix)
  }

