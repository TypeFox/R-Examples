construct.s <-
function (ncoef, lambda, fitting.method, extrapolation = NULL)
{
nl <- length(lambda)
switch(fitting.method,
"quad" = ngamma <- 3,
"line" = ngamma <- 2,
"nonl" = ngamma <- 3,
"logl" = ngamma <- 2,
"log2"= ngamma <- 2
)
# construct a matrix of 0
null.mat <- matrix(0, nrow = ngamma, ncol = ncoef)
s <- list()
for (j in 1:ncoef) {
# define the first matrix
switch(fitting.method,
"quad" = gamma.vec <- c(-1, -lambda[1], -lambda[1]^2),
"line" = gamma.vec <- c(-1, -lambda[1]),
"nonl" = gamma.vec <- c(1,
1 / (coef(extrapolation[[j]])[3] + lambda[1]),
-coef(extrapolation[[j]])[2] /
(coef(extrapolation[[j]])[3] + lambda[1])^2),
"logl" = gamma.vec <- c(exp(coef(extrapolation)[1, j] +
coef(extrapolation)[2, j] * lambda[1]),
exp(coef(extrapolation)[1, j] + coef(extrapolation)[2, j] *
lambda[1]) * lambda[1]),
"log2" = gamma.vec <- c(exp(coef(extrapolation[[j]])[1] +
coef(extrapolation[[j]])[2] * lambda[1]),
exp(coef(extrapolation[[j]])[1] +
coef(extrapolation[[j]])[2] * lambda[1]) * lambda[1])
)
a <- null.mat
# exchange the column j with the deviation of G(gamma, lambda)
a[, j] <- gamma.vec
for (i in 2:nl) {
switch(fitting.method,
"quad" = gamma.vec <- c(-1, -lambda[i], -lambda[i]^2),
"line" = gamma.vec <- c(-1, -lambda[i]),
"nonl" = gamma.vec <- c(1,
1 / (coef(extrapolation[[j]])[3] + lambda[i]),
-coef(extrapolation[[j]])[2] /
(coef(extrapolation[[j]])[3] + lambda[i])^2),
"logl" = gamma.vec <- c(exp(coef(extrapolation)[1, j] +
coef(extrapolation)[2, j] * lambda[i]),
exp(coef(extrapolation)[1, j] + coef(extrapolation)[2, j] *
lambda[i]) * lambda[i]),
"log2" = gamma.vec <- c(exp(coef(extrapolation[[j]])[1] +
coef(extrapolation[[j]])[2] * lambda[i]),
exp(coef(extrapolation[[j]])[1] +
coef(extrapolation[[j]])[2] * lambda[i]) * lambda[i])
)
b <- null.mat
b[, j] <- gamma.vec
a <- cbind(a, b)
}
s[[j]] <- a
}
s <- t(matrix(unlist(lapply(s, t), recursive = FALSE), nrow = nl * ncoef,
ncol = ngamma * ncoef))
return(s)
}

