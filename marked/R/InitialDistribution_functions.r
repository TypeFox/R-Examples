#' HMM Initial state distribution functions
#' 
#' Functions that compute the initial probability distribution for the states. 
#' Currently only CJS, MS models and MS models with state uncertainty are included and
#' these all use cjs_delta to assign a known state.
#'  
#' @param pars list of real parameter matrices (id by occasion) for each type of parameter
#' @param m number of states
#' @param F initial occasion vector 
#' @param T number of occasions
#' @param start matrix with values that are first occasion and for some CJS type models the state of first observation
#' @aliases cjs_delta 
#' @export
#' @return 2-d array of initial state probability vectors for each id
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. 275p. 
cjs_delta=function(pars,m,F,T,start) 
{
	if(is.list(m))m=m$ns*m$na+1
    delta=matrix(0,nrow=nrow(start),ncol=m)
	if(any(is.na(start[,1])))stop("Unknown initial state")
	delta[cbind(1:nrow(start),start[,1])]=1
	delta 
}
