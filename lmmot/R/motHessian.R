##' Hessian matrix of log-Likelihood for right censored Multiple Ordinal Tobit (MOT) model.
##'
##' @title Hessian matrix of log-Likelihood for mot model
##'
##' @param param parameter vector: (beta_0, beta_1, ... , beta_m, sigma).
##' @param xx design matrix of the model.
##' @param y observation vector.
##' @param tau threshold vector from tau_1 to tau_K.
##'
##' @return hessian matrix, summarized over all observations.
##'
##' @seealso \link[lmmot]{lmmot} 
##' @author Marvin Wright

motHessian <- function(param,xx,y,tau) {
  x <- xx
	
	sigma <- param[length(param)]
	beta <- param[-length(param)]
	
	n <- length(y)
	m <- length(param)
	K <- length(tau)
	
	yy <- t(x) %*% beta
	hess <- matrix(0,m,m)
	
	# non censored data
	index <- y < tau[1]
	if (sum(index) > 0) {
	
		# d^2/dbeta^2
		dbetabeta <- -1/sigma^2 * (x[,index] %*% t(x[,index]))
		
		# d^2/dsigma^2
		dsigmasigma <- sum(1/sigma^2 - 3/sigma^4 * (y[index] - yy[index])^2)
		
		# d^2/(dsigma*dbeta)
		dsigmabeta <- -2/sigma^3 * x[,index] %*% as.matrix(y[index] - yy[index])

		# construct hessian matrix
		hess[1:(m-1),1:(m-1)] <- dbetabeta
		hess[m,m] <- dsigmasigma
		hess[m,1:(m-1)] <- dsigmabeta
		hess[1:(m-1),m] <- dsigmabeta
	}
	
	# censored data, categories 1..K-1
	if (K > 1) {
		for (k in 1:(K-1)) {
			index <- (y >= tau[k] & y < tau[k+1])
			if (sum(index) > 0) {
			
				vk <- tau[k] - yy[index]
				vkk <- tau[k+1] - yy[index]
				
				dk <- dnorm(vk/sigma)
				dkk <- dnorm(vkk/sigma) 
				
				pk <- pnorm(vk/sigma)
				pkk <- pnorm(vkk/sigma)
				
				# d^2/dbeta^2
				temp1 <- 1/sigma * (vk*dk - vkk*dkk)/(pkk - pk)
				temp2 <-((dk - dkk)/(pkk - pk))^2	
				dbetabeta <- 1/sigma^2 * x[,index] %*% (t(x[,index]) * (temp1 - temp2))
				
				# d^2/dsigma^2
				temp1 <- ((vk^3/sigma^2 - 2*vk)*dk - (vkk^3/sigma^2 - 2*vkk)*dkk)/(pkk - pk)
				temp2 <- 1/sigma * ((vk*dk - vkk*dkk)/(pkk - pk))^2
				dsigmasigma <- sum(1/sigma^3 * (temp1 - temp2))
				
				# d^2/(dsigma*dbeta)
				temp1 <- ((vk^2/sigma^2 - 1)*dk - (vkk^2/sigma^2 - 1)*dkk)/(pkk - pk)
				temp2 <- 1/sigma * ((vk*dk - vkk*dkk)*(dk - dkk))/(pkk - pk)^2
				dsigmabeta <- 1/sigma^2 * x[,index] %*% as.matrix(temp1 - temp2)
				
				# add to hessian matrix
				hess[1:(m-1),1:(m-1)] <- hess[1:(m-1),1:(m-1)] + dbetabeta
				hess[m,m] <- hess[m,m] + dsigmasigma
				hess[m,1:(m-1)] <- hess[m,1:(m-1)] + dsigmabeta
				hess[1:(m-1),m] <- hess[1:(m-1),m] + dsigmabeta
			}
		}
	}
	
	# last category (K)
	index <- (y >= tau[K])
	if (sum(index) > 0) {
	
		vk <- tau[K] - yy[index]
		dk <- dnorm(vk/sigma)
		pk <- pnorm(vk/sigma)
		
		# d^2/dbeta^2
		temp1 <- 1/sigma * (vk*dk)/(1 - pk)
		temp2 <-(dk/(1 - pk))^2
		dbetabeta <- 1/sigma^2 * x[,index] %*% (t(x[,index]) * (temp1 - temp2))
		
		# d^2/dsigma^2
		temp1 <- ((vk^3/sigma^2 - 2*vk)*dk)/(1 - pk)
		temp2 <- 1/sigma * ((vk*dk)/(1 - pk))^2
		dsigmasigma <- sum(1/sigma^3 * (temp1 - temp2))
		
		# d^2/(dsigma*dbeta)
		temp1 <- ((vk^2/sigma^2 - 1)*dk)/(1 - pk)
		temp2 <- 1/sigma * (vk*dk*dk)/(1 - pk)^2
		dsigmabeta <- 1/sigma^2 * x[,index] %*% as.matrix(temp1 - temp2)
		
		# add to hessian matrix
		hess[1:(m-1),1:(m-1)] <- hess[1:(m-1),1:(m-1)] + dbetabeta
		hess[m,m] <- hess[m,m] + dsigmasigma
		hess[m,1:(m-1)] <- hess[m,1:(m-1)] + dsigmabeta
		hess[1:(m-1),m] <- hess[1:(m-1),m] + dsigmabeta
	}

	return(hess)
}