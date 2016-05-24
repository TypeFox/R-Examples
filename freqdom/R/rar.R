#' Simulate \code{n} observarions multivariate autoregressive time series, i.e.
#' \deqn{ X_t = \sum_{k=0}^p A_k X_{t-k} + \varepsilon_t,}
#' where \eqn{\varepsilon_t} is a \code{d}-dimensional white noise and \eqn{A_k} are \eqn{d \times d} matrices 
#' and \eqn{X_t = 0} for \eqn{t \leq 0}.
#'
#' @title Simulate a multivariate autoregressive time series
#' @param n number of observations to generate
#' @param d dimension of the process
#' @param Psi serie of regression operators (if one matrix is given it is treated as regressor with lag 1)
#' @param first the first element of a series
#' @param noise the noise we want to add
#' @param sd standard deviation of the gaussian noise if the noise wasn't provided
#' @return an AR series of vectors
#' @importFrom graphics plot title

#' @export
#' @examples
#' nbase = 10
#' Psi = t((1:nbase) %*% t(sin(1:nbase * 2*pi/nbase)) / (nbase*nbase))
#' process = rar(30, Psi=Psi, sd=0.2)
#' pdf(file='simulated.arh1.pdf')
#' plot(process)
#' title("Simulated ARH(1)")
#' dev.off()
rar = function(n, d=NULL, Psi=NULL, first=NULL, noise=NULL, sd=1)
{
	if (is.null(d))
		d = dim(Psi)[1]
	if (is.null(d))
		d = 2
  
	if (!is.positiveint(n))
	  stop ("n must be a positive integer")
	if (!is.positiveint(d))
	  stop ("d must be a positive integer")
	
	nbasis = d
	
	if (is.null(nbasis))
		stop("Can't determine the dimension. Specify d or give Psi.")
	if (nbasis < 1)
		stop("Wrong dimension. d must be positive.")
	if (sd < 0)
		stop("Wrong standard deviation. sd must be positive.")
	
	# if no operator Psi is specified then use Identity
	if (is.null(Psi))
		Psi = diag(nbasis)

	# build coefficients matrix (initially null)
	coef = matrix(0,nbasis,n)

	if (is.null(noise)){
		noise = function(nlen) {
			rnorm(nlen, sd=sd)
		}
	}

	# first element is just random 
	if (is.null(first)){
		first = noise(nbasis)
	}

	D = dim(Psi)[3]
	if (is.na(D)){
		Psi = array(Psi, c(1,dim(Psi)))
		D=1
	}

	coef[,1] = first
	# each next follow Y = Psi(X) + noise
	for (i in 2:n){
		last = min(D,i-1)
		for (j in 1:last){
			coef[,i] = coef[,i] + as.matrix(Psi[j,,]) %*% coef[,i-j]
		}
		coef[,i] = coef[,i] + noise(nbasis)
	}

	# return time series
	t(coef)
}
