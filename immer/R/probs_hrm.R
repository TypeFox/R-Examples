

#######################################################
# Probabilities in hierarchical rater model
probs_hrm <- function( x , xi , phi , psi , K , x_ind = NULL , useRcpp=FALSE){
	if ( ! useRcpp ){
		probs <- probs_hrm_R( x , xi , phi , psi , K , x_ind )
					} else {
		probs <- .Call("probs_hrm_rcpp" , x , xi , phi , psi , K , x_ind ,
						PACKAGE="immer")
					}

    return( probs)
        }
####################################################		

# SEXP probs_hrm_rcpp( SEXP x_, SEXP xi_, SEXP phi_, SEXP psi_, SEXP K_, SEXP x_ind_ ){
# BEGIN_RCPP         
#     // # probs_hrm <- function( x , xi , phi , psi , K , x_ind = NULL ){         
#     Rcpp::NumericVector x(x_);          
#     Rcpp::NumericVector xi(xi_);  
#     Rcpp::NumericVector phi(phi_);  
#     Rcpp::NumericVector psi(psi_);  
#     int K = as<int>(K_);  
#     Rcpp::NumericVector x_ind(x_ind_);  

######################################################################
probs_hrm_R <- function( x , xi , phi , psi , K , x_ind = NULL ){
    N <- length(xi)	
	KM <- matrix( 0:K , nrow=N , ncol=K+1 , byrow=TRUE )
    p1 <- exp( - ( KM - ( xi + phi ) )^2 / ( 2 * psi ) )	
    probs <- p1 / rowSums(p1 , na.rm=TRUE)
	if ( ! is.null(x) ){
		ind <- cbind( 1:N , x+1 )
		probs <- probs[ind ]
					}
	if ( ! is.null( x_ind) ){
		probs <- ifelse( x_ind == 0 , 1 , probs )
							}
    return( probs)
        }
####################################################	