wsvm.kernel <- function(X, U, kernel = list(type = 'linear', par = NULL)) {
# Description 
#     Compute kernel K(X, U)
# Usage
#     K = svm.kernel(X, U, kernel)
# Input
#     X = X matrix
#     U = U matrix
#     kernel = kernel list
#        kernel$kind = type of kernel 
#           eg. 'linear', 'poly' and 'rbf'
#          kernel$par = parameter of kernel 
#           eg. par = degree for 'poly' and  par = scale for 'rbf'
# Output
#     K = nrow(X) by nrow(U) kernel matrix
   if(!is.matrix(X)) X <- as.matrix(X)
   if(!is.matrix(U)) U <- as.matrix(U)
   
   if (kernel$type == 'linear') K <- (X %*% t(U))
   else if (kernel$type == 'poly') K <- (1 + X %*% t(U))^kernel$par
   else if (kernel$type == 'rbf'){
      a = as.matrix(apply(X^2, 1, 'sum'))     
      b = as.matrix(apply(U^2, 1, 'sum'))
      one.a = matrix(1, ncol = nrow(b))     
      one.b = matrix(1, ncol = nrow(a))
      K1 = one.a %x% a
      K2 = X %*% t(U)
      K3 = t(one.b %x% b)
      K = exp(-(K1 - 2 * K2 + K3)/(2 * kernel$par^2))
   }
   return(K)
}     
