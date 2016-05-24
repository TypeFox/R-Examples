logRegSubspace <- function(D, Y, A, JJs, q){

mle <- matrix(glm.fit(Y, D, family = binomial(link = "logit"))$coefficients, ncol = 1)

## in case of no violations of constraints just compute
## unconstrained estimates
if (length(A) == 0){beta <- mle} 

## check whether we only have ordered factors and
## all constraints are active. In this case,
## beta is then simply a vector of 0's
all.act <- identical(unlist(JJs), A)

if ((length(A) > 0) && (all.act == 1)){
    beta <- matrix(0, nrow = length(mle))
} # end if

# ================================
if ((length(A) > 0) && (all.act == 0)){

res <- shrinkBeta(Y, A, JJs, q)
Y.col <- res$Y.col
sums <- res$sums
rems <- res$rems
JJs.A <- res$JJs.A

## compute unconstrained estimator on subspace
beta.col <- matrix(glm.fit(Y.col, D, family = binomial(link = "logit"))$coefficients, ncol = 1)

## expand estimate to receive estimate in original dimension back
beta <- expandBeta(beta.col, sums, JJs.A)$beta

} # end if

return(list("beta" = beta))
}













#
