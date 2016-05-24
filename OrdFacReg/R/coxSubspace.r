coxSubspace <- function(ttf, tf, Y, A, JJs, q){

mle <- matrix(eha::coxreg.fit(Y, Surv(ttf, tf), max.survs = length(tf), strats = rep(1, length(tf)))$coefficients, ncol = 1)

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

    mini <- min(unlist(JJs))
    res <- shrinkBeta(Y, A, JJs, q)
    Y.col <- res$Y.col
    sums <- res$sums
    rems <- res$rems
    JJs.A <- res$JJs.A

    ## compute unconstrained estimator on subspace
    beta.col <- matrix(eha::coxreg.fit(Y.col, Surv(ttf, tf), max.survs = length(tf), strats = rep(1, length(tf)))$coefficients, ncol = 1)

    ## expand estimate to receive estimate in original dimension back
    beta <- expandBeta(beta.col, sums, JJs.A)$beta

} # end if

return(list("beta" = beta))
}













#
