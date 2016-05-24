##' Fisher information matrix for right censored Multiple Ordinal Tobit (MOT) model.
##'
##' @title Fisher information for mot model
##'
##' @param param parameter vector: (beta_0, beta_1, ... , beta_m, sigma).
##' @param xx design matrix of the model.
##' @param tau threshold vector from tau_1 to tau_K.
##'
##' @return fisher information matrix, summarized over all observations.
##' 
##' @seealso \link[lmmot]{lmmot} 
##' @author Marvin Wright

motFisher <- function(param,xx,tau) {
  x <- xx
  
  sigma <- param[length(param)]
  beta <- param[-length(param)]
  
  m <- length(param)
  K <- length(tau)
  
  yy <- t(x) %*% beta
  fish <- matrix(0,m,m)
  
  # non censored data
  p <- pnorm((tau[1]-yy)/sigma)

  # d^2/dbeta^2
  dbetabeta <- 1/sigma^2 * x %*% (t(x) * as.vector(p))
    
  # d^2/dsigma^2
  dsigmasigma <- sum(2/sigma^2 * as.vector(p))
  
  # d^2/(dsigma*dbeta)
  dsigmabeta <- 0
  
  # construct fisher info
  fish[1:(m-1),1:(m-1)] <- dbetabeta
  fish[m,m] <- dsigmasigma
  fish[m,1:(m-1)] <- dsigmabeta
  fish[1:(m-1),m] <- dsigmabeta
  
  # censored data, categories 1..K-1
  if (K > 1) {
    for (k in 1:(K-1)) {
      
      vk <- tau[k] - yy
      vkk <- tau[k+1] - yy
      
      dk <- dnorm(vk/sigma)
      dkk <- dnorm(vkk/sigma) 
      
      pk <- pnorm(vk/sigma)
      pkk <- pnorm(vkk/sigma)
      
      # d^2/dbeta^2
      temp1 <- 1/sigma * (vk*dk - vkk*dkk)
      temp2 <- (dk - dkk)^2/(pkk - pk) 
      temp <- as.vector(temp1 - temp2)
      temp[pkk-pk == 0] <- 0
      dbetabeta <- -1/sigma^2 * x %*% (t(x) * temp)
      
      # d^2/dsigma^2
      temp1 <- ((vk^3/sigma^2 - 2*vk)*dk - (vkk^3/sigma^2 - 2*vkk)*dkk)
      temp2 <- 1/sigma * (vk*dk - vkk*dkk)^2/(pkk - pk)
      temp <- as.vector(temp1 - temp2)
      temp[pkk-pk == 0] <- 0
      dsigmasigma <- -sum(1/sigma^3 * temp)
      
      # d^2/(dsigma*dbeta)
      temp1 <- (vk^2/sigma^2 - 1)*dk - (vkk^2/sigma^2 - 1)*dkk
      temp2 <- 1/sigma * ((vk*dk - vkk*dkk)*(dk - dkk))/(pkk - pk)
      temp <- as.vector(temp1 - temp2)
      temp[pkk-pk == 0] <- 0
      dsigmabeta <- -1/sigma^2 * x %*% temp
      
      # add to fisher info
      fish[1:(m-1),1:(m-1)] <- fish[1:(m-1),1:(m-1)] + dbetabeta
      fish[m,m] <- fish[m,m] + dsigmasigma
      fish[m,1:(m-1)] <- fish[m,1:(m-1)] + dsigmabeta
      fish[1:(m-1),m] <- fish[1:(m-1),m] + dsigmabeta
      
    }
  }
  
  # last category (K)
  vk <- tau[K] - yy
  dk <- dnorm(vk/sigma)
  pk <- pnorm(vk/sigma)
  
  # d^2/dbeta^2
  temp1 <- 1/sigma * (vk*dk)
  temp2 <- (dk^2/(1 - pk))
  temp <- as.vector(temp1 - temp2)
  temp[pk == 1] <- 0
  dbetabeta <- -1/sigma^2 * x %*% (t(x) * temp)
  
  # d^2/dsigma^2
  temp1 <- ((vk^3/sigma^2 - 2*vk)*dk)
  temp2 <- 1/sigma * ((vk*dk)^2/(1 - pk))
  temp <- as.vector(temp1 - temp2)
  temp[pk == 1] <- 0
  dsigmasigma <- -sum(1/sigma^3 * temp)
  
  # d^2/(dsigma*dbeta)
  temp1 <- (vk^2/sigma^2 - 1)*dk
  temp2 <- 1/sigma * (vk*dk*dk)/(1 - pk)
  temp <- as.vector(temp1 - temp2)
  temp[pk == 1] <- 0
  dsigmabeta <- -1/sigma^2 * x %*% temp
  
  # add to fisher info
  fish[1:(m-1),1:(m-1)] <- fish[1:(m-1),1:(m-1)] + dbetabeta
  fish[m,m] <- fish[m,m] + dsigmasigma
  fish[m,1:(m-1)] <- fish[m,1:(m-1)] + dsigmabeta
  fish[1:(m-1),m] <- fish[1:(m-1),m] + dsigmabeta
    
  return(fish)
}