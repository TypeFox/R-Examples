Mimpute <- function(data, patused, mu, sig) {
# This function imputes the missing data, given the observed and mu and sig
# for a single pattern it goeas through each patterns and uses the
# conditional distribution of missing given observed and mu and sig to
# impute from the appropriate posterior distribution 
 ni <- nrow(data)
 pp <- ncol(data)
 indm <- which(is.na(patused))
 indo <- which(!is.na(patused))
 pm <- length(indm)
 po <- length(indo)
 muo <- mu[indo]
 mum <- mu[indm]
 sigooi <- solve(sig[indo, indo])
 sigmo <- matrix(sig[indm, indo], pm, po)
 sigmm <- matrix(sig[indm, indm], pm, pm)
 ss1 <- sigmo %*% sigooi
 varymiss <- sigmm - ss1 %*% t(sigmo)
 expymiss <- matrix(mum, ni, pm, byrow = TRUE) + 
             (data[, indo] - matrix(muo, ni, po, byrow = TRUE)) %*% t(ss1)
 if (pm == 1) {
    a <- sqrt(varymiss)
 } else {
    svdvar <- svd(varymiss)
    a <- diag(sqrt(svdvar$d)) %*% t(svdvar$u)
 }
 data[, indm]= matrix(rnorm(ni * pm), ni, pm) %*% a + expymiss
 data
}
