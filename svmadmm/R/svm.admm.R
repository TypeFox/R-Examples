#' linear/nonlinear svm solver based ADMM and IADMM algorithms
#' 
#' @description
#'   \code{svm.admm} is a simple function for solving large-scale regularized linear/nonlinear
#'   classification by using ADMM and IADMM algorithms. This function provides 
#'   linear L2-regularized primal classification (both ADMM and IADMM are available), 
#'   kernel L2-regularized dual classification (IADMM) as well as L1-regularized primal
#'   classification (both ADMM and IADMM are available). The training of the models perform well 
#'   practice.
#'
#' @param x.tr a n*p data matrix. Each row stands for an example (sample, point) and each column stands for a dimension (feature, variable).
#' @param y.tr a n-length vector. The values correspond to class labels. 
#' @param type \code{svmadmm} can provide 3 types of linear/kernel models. 
#'             The default value for \code{type} is 0. See details below. Valid options are:
#'   \itemize{  
#'         \item 0 -- L2-regularized kernel svm (dual)
#'         \item 1 -- L2-regularized linear svm (primal)
#'         \item 2 -- L1-regularized linear svm (primal)
#'     }
#' @param kernel the kernel used in training and predicting when \code{type} is 0.
#'             The default value for \code{kernel} is radial. See details below. Valid options are:
#'   \itemize{  
#'         \item \code{radial} -- The Gaussian RBF kernel k(x,x') = exp(-sigma \|x - x'\|^2) 
#'         \item \code{linear} -- The Linear kernel k(x,x') = <x, x'>
#'         \item \code{polynomial} -- The Polynomial kernel k(x,x') = (scale <x, x'> + offset)^{degree}
#'     }
#' @param sigma The inverse kernel width used by the Gaussian.
#' @param degree The degree of the polynomial kernel function. This has to be an positive integer.
#' @param scale The scaling parameter of the polynomial kernel is a convenient way of normalizing patterns without the need to modify the data itself.
#' @param offset The offset used in a polynomial kernel.
#' @param algo the algorithm to solve the problem w.r.t. \code{type}.
#' @param lambda regularization constant (default: 1). Rules the trade-off
#'   between regularization and correct classification on \code{data}.
#' @param rho regularization constant (default: 1).
#' @param eps epsilon in the termination condition.
#' @details
#' \code{svmadmm} internally computing kernel matrix when \code{type} is 0, which is based by the package \strong{kernlab}.
#' @return
#' An list containing the fitted model, including:
#' \item{alpha}{ A solution for dual form svm classification.}
#' \item{beta}{ A solution for primal form svm classifciation, also are the model weights.}
#' \item{type}{ An integer correspinding to \code{type}.}
#' \item{kernel}{ A function to define the kernel.}
#' \item{x.tr}{ The training input data.}
#' \item{y.tr}{ The training output data.}
#' @examples
#' library(svmadmm)
#' n = 100
#' p = 10
#' x = matrix(runif(2 * n * p, -1, 1), nrow = 2 * n)
#' y = sign(x[, 1])
#' y.ind = sample(1 : (2 * n), n / 10, replace = FALSE)
#' y[y.ind] = - y[y.ind]
#' x.tr = x[1 : n, ]
#' y.tr = y[1 : n]
#' x.te = x[-(1 : n), ]
#' y.te = y[-(1 : n)]
#' model = svm.admm(x.tr, y.tr)
#' fit = svm.predict(x.te, model)

#' @export 
svm.admm <- function(x.tr, y.tr, type = 0, kernel = "radial", sigma = 1 / ncol(x.tr), degree = 1, scale = 1, offset = 1, algo = "iadmm", lambda = 1, rho = 1, eps = 1e-2)
{
	# if (!requireNamespace("kernlab")) {
	# 	# warning(The package kernlab was not installed)
	# 	# install.packages("kernlab")
 #  		stop('The package kernlab was not installed')
	# }
	# L2-regularized kernel svm (dual)
	if(type == 0)
	{
		if (algo == "admm"){
			warning("admm can't solve the problem efficiently, the solution is solved by iadmm")
		} 

		if (kernel == "radial")
		{
			rbf1 = kernlab::rbfdot(sigma = sigma)
		} 
		else if (kernel == "linear")
		{
			rbf1 = kernlab::vanilladot()

		} 
		else if (kernel == "polynomial")
		{
			rbf1 = kernlab::polydot(degree = degree, scale = scale, offset = offset)
		}
		  n = nrow(x.tr)
		  p = ncol(x.tr)
		  k1 = kernlab::kernelPol(rbf1, x.tr, z = y.tr)

		  rbf2 = kernlab::vanilladot()
		  x = matrix(rep(1, n * 1), n, 1)
		  k2 = kernlab::kernelPol(rbf2, x, z = y.tr)
		  k = k1 + k2

		  zeta = max(apply(abs(k), 2, sum))
		  iadmm <- list(alpha = rep(0, n), eta = rep(0, n), mu = rep(0, n))
          old <- list(alpha = rep(1, n), eta = rep(1, n), mu = rep(1, n))
          r = 100
		  s = 100
		  eps = eps
		  # e.r = eps
		  # e.s = eps
          c = 1 / lambda
          ratio = 1
          # com3 = norm( as.matrix(rep(1, n)), type = "f")

          
		while ( ratio > eps ){
		    iadmm$eta = pmax(c - iadmm$alpha - iadmm$mu, 0)
		    iadmm$mu = old$mu + ( - c + iadmm$alpha + iadmm$eta)
		    iadmm$alpha = as.vector( 1 / (zeta + rho) * (zeta * iadmm$alpha - k %*% iadmm$alpha + 1 - rho * (iadmm$eta - c + iadmm$mu)))
		    iadmm$alpha = pmax(iadmm$alpha, 0)
		    ratio = sum(abs(iadmm$alpha - old$alpha)) / sum(abs(iadmm$alpha))
		    old = iadmm
		  
			# r = norm( as.matrix(iadmm$alpha + iadmm$eta - 1), type = "f")
			# s = rho * norm( as.matrix( iadmm$eta - old$eta ), type = "f")
			# com1 = norm( as.matrix(iadmm$alpha), type = "f")
			# com2 = norm( as.matrix(iadmm$eta), type = "f")
			# e.r = sqrt(n) * eps + eps * max(c(com1, com2, com3))
			# e.s = sqrt(p) * eps + eps * rho * norm( as.matrix(iadmm$mu), type = "f")
		}

		  list(alpha = iadmm$alpha, type = type, kern = rbf1, x.tr = x.tr, y.tr = y.tr)
	
		}


    # L2-regularized linear svm (primal)
	else if (type == 1)
	{ 
		q = 2
		if (algo == "admm"){
		  n = nrow(x.tr)
		  x.tr = cbind(rep(1, n), x.tr)
		  p = ncol(x.tr)
		  tyx = t(y.tr * x.tr)
		  inv = solve(tyx %*% t(tyx) + diag(2 * lambda / rho, p))
		  output = .C("linearadmm", 
		  	        as.integer(n), 
		  	        as.integer(p), 
		  	        as.double(x.tr), 
		  	        as.double(y.tr), 
		  	        as.double(tyx), 
		  	        as.double(inv), 
		  	        as.double(lambda), 
		  	        as.double(rho), 
		  	        as.integer(q),
		  	        as.double(eps),
		  	        output = double(p),		  	        
		  	        PACKAGE = "svmadmm"
		  	        )
		  beta = output$output
		  list(beta = beta, type = type, x.tr = x.tr, y.tr = y.tr)
		} 

		else if (algo == "iadmm"){		
		  n = nrow(x.tr)
		  x.tr = cbind(rep(1, n), x.tr)
		  p = ncol(x.tr)
		  tyx = t(y.tr * x.tr)
		  zeta = rho * sum(diag(t(tyx) %*% (tyx)))
		  output = .C("lineariadmm", 
		  			as.integer(n), 
		  			as.integer(p), 
		  			as.double(x.tr), 
		  			as.double(y.tr), 
		  			as.double(tyx), 
		  			as.double(zeta), 
		  			as.double(lambda), 
		  			as.double(rho), 
		  			as.integer(q),
		  			as.double(eps),
		  			output = double(p),
		  			PACKAGE = "svmadmm")
		  beta = output $ output
		  list(beta = beta, type = type, x.tr = x.tr, y.tr = y.tr)
		}
		else{
			stop("parameter algo can only be admm or iadmm.")
		}
    
	} 
	# L1-regularized linear svm (primal)
	else if (type == 2)
	{
		q = 1
		if (algo == "admm"){
		  warning("admm don't really solve the problem efficiently, the solution is solved by iadmm")
		}

		  n = nrow(x.tr)
		  x.tr = cbind(rep(1, n), x.tr)
		  p = ncol(x.tr)
		  tyx = t(y.tr * x.tr)
		  zeta = rho * sum(diag(t(tyx) %*% (tyx)))
		  output = .C("lineariadmm", 
		  			as.integer(n), 
		  			as.integer(p), 
		  			as.double(x.tr), 
		  			as.double(y.tr), 
		  			as.double(tyx), 
		  			as.double(zeta), 
		  			as.double(lambda), 
		  			as.double(rho), 
		  			as.integer(q),
		  			as.double(eps),
		  			output = double(p),
		  			PACKAGE = "svmadmm")
		  beta = output $ output
		  list(beta = beta, type = type, x.tr = x.tr, y.tr = y.tr)
	}
}