#
#

sar <- function(X, y, W, burnin=5000, Nsim=10000) {
    detval <- lndet_imrw(W)
    result<- .Call('HSAR_sar_cpp_arma', PACKAGE = 'HSAR', X, y, W, detval, burnin, Nsim)
    class(result) <- "mcmc_sar"
    result
}
