#### from devtools::use_rcpp()
#' @useDynLib chopthin
#' @importFrom Rcpp sourceCpp

a <- function(x){x} ## dummy to ensure that there is SOME R code;
                    ## seems to be necessary to ensure that the above
                    ## commands to import the dynamic libraries are
                    ## actually being processed.
