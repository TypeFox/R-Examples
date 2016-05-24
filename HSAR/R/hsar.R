#
#

hsar <- function(X, y, W=NULL, M=NULL, Z, Unum, burnin=5000, Nsim=10000) {
    if( is.null(W) & is.null(M) ){
      cat("Both weight matrices can not be NULL")
      return(NULL)
    }
    # Special case where rho =0 ; dependent regional effect 
    if (is.null(W)){
      detval <- lndet_imrw(M)
      result <- .Call('HSAR_hsar_cpp_arma_rho_0', PACKAGE = 'HSAR', X, y, M, Z, detval, Unum, burnin, Nsim)
      class(result) <- "mcmc_hsar_rho_0"
    }
    # Special case where lamda =0 ; independent regional effect
    if ( is.null(M)){
      detval <- lndet_imrw(W)
      result <- .Call('HSAR_hsar_cpp_arma_lambda_0', PACKAGE = 'HSAR', X, y, W, Z, detval, Unum, burnin, Nsim)
      class(result) <- "mcmc_hsar_lambda_0"
    }
    # Full HSAR model
    if ( (!is.null(M)) & (!is.null(W))){
      detval <- lndet_imrw(W)
      detvalM <- lndet_imrw(M)
      result <- .Call('HSAR_hsar_cpp_arma', PACKAGE = 'HSAR', X, y, W, M, Z, detval, detvalM, Unum, burnin, Nsim)
      class(result) <- "mcmc_hsar"
    }
    
    return(result)
}
