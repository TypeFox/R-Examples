##' lgcppars function
##'
##' A function for setting the parameters sigma, phi and theta for \code{lgcpPredict}. Note that the returned
##' set of parameters also features mu=-0.5*sigma^2, gives mean(exp(Y)) = 1.
##'
##' @param sigma sigma parameter
##' @param phi phi parameter
##' @param theta this is 'beta' parameter in Brix and Diggle (2001)
##' @param mu the mean of the latent field, if equal to NULL, this is set to -sigma^2/2
##' @param beta ONLY USED IN case where there is covariate information.
##' @seealso \link{lgcpPredict}
##' @export

lgcppars <- function(sigma=NULL,phi=NULL,theta=NULL,mu=NULL,beta=NULL){
    
    if(is.null(mu)){
	    return(list(sigma=sigma,phi=phi,mu=-0.5*sigma^2,theta=theta,beta=beta))
	}
	else{
	    return(list(sigma=sigma,phi=phi,mu=mu,theta=theta,beta=beta))
	}	
	# NB choice of parameter mu=-0.5*sigma^2 gives mean(exp(Y)) = 1  					
}						
