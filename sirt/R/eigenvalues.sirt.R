

################################################
# calculation of first D eigenvalues
eigenvalues.sirt <- function (X,D, maxit=200 , conv=10^(-6) ){ 
	.Call("eigenvaluesDsirt", X,D,maxit,conv, PACKAGE = "sirt")
					}				
# extern "C" {
# SEXP eigenvalues_sirt( SEXP X_, SEXP D_, SEXP maxit_, SEXP conv_) ;
# }
