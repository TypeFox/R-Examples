LSEsubspace <- function(D, Y, A, JJs, q){

lse <- lmLSE(D, Y)$beta

## in case of no violations of constraints just compute
## unconstrained estimates
if (length(A) == 0){beta <- lse} 

## check whether we only have ordered factors and
## all constraints are active. In this case,
## beta is then simply a vector of 0's
all.act <- identical(unlist(JJs), A)

if ((length(A) > 0) && (all.act == 1)){
    beta <- matrix(0, nrow = length(lse))
} # end if

# ================================
if ((length(A) > 0) && (all.act == 0)){

res <- shrinkBeta(Y, A, JJs, q)
Y.col <- res$Y.col
sums <- res$sums
rems <- res$rems
JJs.A <- res$JJs.A


## compute unconstrained estimator on subspace
mat <- t(Y.col) %*% Y.col
beta.col <- MASS::ginv(mat) %*% t(Y.col) %*% D


## expand estimate to receive estimate in original dimension back
beta <- expandBeta(beta.col, sums, JJs.A)$beta

} # end if

return(list("beta" = beta))
}













#
