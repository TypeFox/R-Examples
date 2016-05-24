#' Computation of the Empiric Cramer'v.
#' 
#' This function computes the Empiric Cramer's V for a binary data set.
#' @param x a binary matrix.
#' @return Return the matrix of the Empiric Cramer's V.
#' 
#' @export
#' 
ComputeEmpiricCramer <- function(x){
  if ( (is.matrix(x)==FALSE) || any((x==0) + (x==1) == 0) )
    stop("The input parameter x must be a binary matrix")
  
  alpha <- colMeans(x)
  VcramerEmpiric <- matrix(0, ncol(x), ncol(x))
  rownames(VcramerEmpiric) <- colnames(VcramerEmpiric) <- colnames(x)
  for (j in 1:ncol(x)){
    obs <- rbind((t(x[,j])%*%(x)), (t(x[,j])%*%(1-x)), (t(1-x[,j])%*%(x)), (t(1-x[,j])%*%(1-x)))/nrow(x)
    th <- rbind(alpha[j]*alpha, alpha[j] * (1-alpha), (1-alpha[j])*alpha, (1-alpha[j])*(1-alpha))
    VcramerEmpiric[j,] <- sqrt(colSums((obs - th)**2 / th))
  }
  return(VcramerEmpiric)
}

#' Computation of the model Cramer'v.
#' 
#' This function computes the model Cramer's V for a binary data set.
#' @param results an instance of S4 class MvBinaryResult (provided by the function MvBinaryEstim)
#' @return Return the matrix of the Empiric Cramer's V.
#' 
#' @export
#' 
ComputeMvBinaryCramer <- function(results){
  if (class(results)!= "MvBinaryResult")
    stop("The imput parameter must be an instance of S4 class MvBinaryResult (provided by the function MvBinaryEstim) ")
  Vcramer <- matrix(0, length(results@blocks), length(results@blocks))
  diag(Vcramer) <- 1
  rownames(Vcramer) <- colnames(Vcramer) <- names(results@blocks)
  for (b in 1:max(results@blocks)){
    who <- which(results@blocks==b)
    for (j in who){
      for (j2 in who){
        if (j!=j2){
          alpha <- c(results@alpha[j], results@alpha[j2])
          epsilon <- c(results@epsilon[j], results@epsilon[j2])
          delta <- c(results@delta[j], results@delta[j2])
          beta <- alpha*delta + (1-alpha)*(1-delta)
          lambda <- (1-epsilon)* alpha + epsilon*delta
          mu <- (1-epsilon)* alpha + epsilon*(1-delta)
          ord <- order(beta)
          probamodel <- rep(0, 4)
          probamodel[1] <- beta[ord[1]] * lambda[ord[1]] * lambda[ord[2]] + (beta[ord[2]]-beta[ord[1]]) * mu[ord[1]] * lambda[ord[2]] + (1- beta[ord[2]]) * mu[ord[1]] * mu[ord[2]]
          probamodel[2] <- beta[ord[1]] * (1-lambda[ord[1]]) * lambda[ord[2]] + (beta[ord[2]]-beta[ord[1]]) * (1-mu[ord[1]]) * lambda[ord[2]] + (1- beta[ord[2]]) * (1-mu[ord[1]]) * mu[ord[2]]
          probamodel[3] <- beta[ord[1]] * lambda[ord[1]] * (1-lambda[ord[2]]) + (beta[ord[2]]-beta[ord[1]]) * mu[ord[1]] * (1-lambda[ord[2]]) + (1- beta[ord[2]]) * mu[ord[1]] * (1-mu[ord[2]])
          probamodel[4] <- beta[ord[1]] * (1-lambda[ord[1]]) * (1-lambda[ord[2]]) + (beta[ord[2]]-beta[ord[1]]) * (1-mu[ord[1]]) * (1-lambda[ord[2]]) + (1- beta[ord[2]]) *(1-mu[ord[1]]) * (1-mu[ord[2]])
          probath <- rep(0, 4)
          probath[1] <- alpha[ord[1]] * alpha[ord[2]]
          probath[2] <- (1-alpha[ord[1]]) * alpha[ord[2]]
          probath[3] <- alpha[ord[1]] * (1-alpha[ord[2]])
          probath[4] <- (1-alpha[ord[1]]) * (1-alpha[ord[2]])
          Vcramer[j,j2] <- Vcramer[j2,j] <-  sqrt(sum(((probamodel - probath)**2)/probath))
        }
      }
    }
  }
  return(Vcramer)
}