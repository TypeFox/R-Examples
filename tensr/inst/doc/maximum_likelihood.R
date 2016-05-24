## ------------------------------------------------------------------------
library(tensr)
p <- c(10, 10, 10, 10)
X <- array(rnorm(prod(p)),dim = p)

cov_Y <- start_ident(p)
cov_Y[[4]] <- 0.9^abs(outer(1:p[4],1:p[4],"-"))
cov_Y[[4]] <- cov_Y[[4]] / det(cov_Y[[4]])^(1/p[4]) #scale to have unit determinant.
Y <- atrans(array(rnorm(prod(p)),dim = p), lapply(cov_Y, mhalf))


cov_Z <- start_ident(p)
diag(cov_Z[[1]]) <- 1:p[1] / prod(1:p[1])^(1/p[1])
diag(cov_Z[[2]]) <- p[2]:1 / prod(1:p[2])^(1/p[2])
cov_Z[[3]] <- 0.9^abs(outer(1:p[3],1:p[3],"-"))
cov_Z[[3]] <- cov_Z[[3]] / det(cov_Z[[3]])^(1/p[3])
Z <- atrans(array(rnorm(prod(p)),dim = p), lapply(cov_Z, mhalf))

## ------------------------------------------------------------------------
holq_X <- holq(X,mode_rep = 1:4, print_diff = FALSE)
holq_Y <- holq(Y,mode_rep = 1:3, print_diff = FALSE)
holq_Z <- holq(Z,mode_rep = 4, mode_diag = c(1,2), print_diff = FALSE)

mle_X <- mle_from_holq(holq_X)
mle_Y <- mle_from_holq(holq_Y)
mle_Z <- mle_from_holq(holq_Z)

## ------------------------------------------------------------------------
cat("Estimates of scale:\n",
    "From X:", mle_X$sig_mle,"\n",
    "From Y:", mle_Y$sig_mle,"\n",
    "From Z:", mle_Z$sig_mle,"\n"
)

## ------------------------------------------------------------------------
cat("     True:", paste("(", paste(format(diag(cov_Z[[1]]), digits = 2),
                                   collapse = ", "), ")",sep = ""), "\n",
    "Estimate:", paste("(", paste(format(diag(mle_Z[[1]][[1]]), digits = 2),
                                  collapse = ", "), ")", sep = ""), "\n")

## ---- results = "hide"---------------------------------------------------
null_distribution <- lrt_null_dist_dim_same(p = p, null_ident = 1:4,
                                            alt_ident = 1:3, itermax = 500)

sig_k <- holq(X,mode_rep = 1:3, print_diff = FALSE)$sig
sig_h <- holq(X,mode_rep = 1:4, print_diff = FALSE)$sig
lrt_stat_val_X <- lrt_stat(sig_null = sig_h,sig_alt = sig_k,p = p)

sig_k <- holq(Y,mode_rep = 1:3, print_diff = FALSE)$sig
sig_h <- holq(Y,mode_rep = 1:4, print_diff = FALSE)$sig
lrt_stat_val_Y <- lrt_stat(sig_null = sig_h,sig_alt = sig_k,p = p)

## ------------------------------------------------------------------------
p_value_x <- mean(null_distribution > lrt_stat_val_X)
p_value_y <- mean(null_distribution > lrt_stat_val_Y)

cat(" p-value using X:", p_value_x,"\n",
    "p-value using Y:", p_value_y,"\n")

## ------------------------------------------------------------------------
par(mar = c(2.5, 2.5, 2, 0), mgp = c(1.5, .5, 0))
hist(null_distribution, xlab = "LRT Stat", main = "Null Distribution")
abline(v = lrt_stat_val_X, col = 2, lwd = 2, lty = 2)
abline(v = lrt_stat_val_Y, col = 4, lwd = 2, lty = 2)
legend("topright", c("From X","From Y"), col = c(2,4), lty = 2, lwd = 2, cex = 0.6)

## ---- results = "hide"---------------------------------------------------
null_distribution <- lrt_null_dist_dim_same(p = p,null_ident = c(2,4),
                                            alt_ident = 4, null_diag = 1,
                                            alt_diag = c(1,2), itermax = 500)

sig_k <- holq(Z, mode_rep = 4, mode_diag = c(1,2), print_diff = FALSE)$sig
sig_h <- holq(Z, mode_rep = c(2,4), mode_diag = 1, print_diff = FALSE)$sig
lrt_stat_val_Z <- lrt_stat(sig_null = sig_h, sig_alt = sig_k, p = p)

p_value_z <- mean(null_distribution > lrt_stat_val_Z)

## ------------------------------------------------------------------------
cat("P-value:", p_value_z,"\n")

