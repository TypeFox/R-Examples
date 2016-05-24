###################################################################################################
# SVM (Support Vector Machine) Version 1.6
# SVM consists of the following functions:
#     svm: main program of SVM
#      svm.compact: fit an SVM model in a compact form
#      svm.predict: predict new.Y at new.X and compute error rate
#      find.nonzero: find Amat.compact and Aind for QP.solve.compact
#     kmeans.process: use K-means for clustering and process response
###################################################################################################
wsvm <- function(X, Y, c.n, kernel = list(type = 'linear', par = NULL), C = 1, eps = 1e-10){
# Description 
#     Fit an svm model
# Usage
#     S = wsvm(X, Y, c.n, kernel, C, eps)
# Input
#     X = input variable matrix
#     Y = output variable vector which will be declared as a matrix in SVM
#     c.n = parameter weight term
#     kernel = kernel list
#        kernel$kind = type of kernel 
#           eg. 'linear', 'poly' and 'rbf'
#          kernel$par = parameter of kernel 
#           eg. par = degree for 'poly' and  par = scale for 'rbf'
#     C = regularization parameter
#     eps = small number
# Output
#     model = a list consists of fit, alpha, bias and sv
#        model$fit = predicted values (n by 1)
#        model$alpha = estimated coefficients (n by 1)
#        model$bias = bias term
#        model$sv = index of support vectors
    argmin <- function(z) {
       index <- 1:length(z)
       argmin <- min(index[z == min(z)])
       return(argmin)
    }
    
    argmax <- function(z){
       index <- 1:length(z)
       argmax <- max(index[z == max(z)])
       return(argmax)
    }
    
   # declare preliminary quantities
   eps <- C * eps
   if(!is.matrix(X)) X <- as.matrix(X)
   if(!is.matrix(Y)) Y <- as.matrix(Y)
   n.data <- nrow(Y) 
   I.n <- diag(rep(1, n.data))   
   Y.n <- Y * c.n   ### Y hat of weighted SVM

   # compute kernel matrix

   K <- wsvm.kernel(X, X, kernel)

   # Solve QP
   # prepare QP
        find.nonzero <- function(Amat){
        # Description 
        #     find Amat.compact and Aind for QP.solve.compact
        # Usage
        #     FN = find.nonzero(Amat)
        # Input
        #     Amat = Amat in solve.QP
        # Output
        #     Amat.compact = compact form of Amat
        #     Aind = indicator of nonzero elements of Amat
           nr <- nrow(Amat)
           nc <- ncol(Amat)
           Amat.compact <- matrix(0, nr, nc)
           Aind <- matrix(0, nr+1, nc)
           for (j in 1:nc){
              index <- (1:nr)[Amat[, j] != 0]
              number <- length(index)
              Amat.compact[1:number, j] <- Amat[index, j]
              Aind[1, j] <- number
              Aind[2:(number+1), j] <- index
           }
           max.number <- max(Aind[1, ])
           Amat.compact <- Amat.compact[1:max.number, ]
           Aind <- Aind[1:(max.number+1), ]
           compact <- list(Amat.compact = Amat.compact, Aind = Aind)
           return(compact)
        }
   Dmat <- K * (Y.n %*% t(Y.n))
   diag(Dmat) <- diag(Dmat) + eps 
   dvec <- c.n             
   Amat <- cbind(Y.n, I.n, -I.n)
   nonzero <- find.nonzero(Amat)
   Amat <- nonzero$Amat.compact
   Aind <- nonzero$Aind  
   bvec <- c(0, rep(0, n.data), rep(-C, n.data))
   # find alpha by QP
   alpha <- solve.QP.compact(Dmat, dvec, Amat, Aind, bvec, meq = 1)$solution

   # compute the index and the number of support vectors 
   S <- 1:n.data
   sv.index <- S[alpha > eps]
   sv.number <- length(sv.index)
   sv.index.C <- S[alpha > eps & alpha < (C - eps)]                 ### &은 and의 의미이고 |는 또는or의 의미이다.
   sv <- list(index = sv.index, number = sv.number, index.C = sv.index.C)
#   print('Number of Support Vectors =')
#   print(sv.number)

   # let alpha = 0 if it is too small
   alpha[-sv.index] <- 0

   # compute bias
   if(length(sv.index.C) == 0){
      sv.index.C <- S[(alpha > eps) & (alpha <= C )]
   }
   bias <- mean(Y[sv.index.C] - K[sv.index.C, sv.index] %*% (alpha[sv.index] * Y[sv.index]))

   # compute fit
   fit <- bias + K %*% (Y * alpha)  
   w <- rep(0,ncol(X))   
   for (i in (1:n.data)) w <- w+alpha[i]*Y[i]*X[i,]

   # prepare output
   model <- list(alpha = alpha, bias = bias, sv = sv, fit = fit, kernel = kernel, w=w, C = C)
   return(model)
}
