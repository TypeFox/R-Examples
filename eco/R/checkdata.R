checkdata <- function(X,Y, supplement, ndim) {
   # check and reorganize inputs 
   if (any(X<0) || any(X>1) || any(Y<0) || any(Y>1))
     stop("Values of X and Y have to be between 0 and 1.")
   ind <- 1:length(X)
   res <- list()
   res$X1type <- res$samp.X1 <- res$X1.W1 <- 0 
   res$X0type <- res$samp.X0 <- res$X0.W2 <- 0  
   
   ## X = 1
   X1.ind <- ind[along=(X==1)]
   if (length(X[X!=1])<length(X)){
      res$X1type <- 1
      res$samp.X1 <- length(X1.ind)
      res$X1.W1 <- Y[X1.ind]
   }

   ## X = 0 
   X0.ind <- ind[along=(X==0)]
   if (length(X[X!=0])<length(X)){
      res$X0type <- 1
      res$samp.X0 <- length(X0.ind)
      res$X0.W2 <- Y[X0.ind]
   }

   XX.ind <- setdiff(ind, union(X0.ind, X1.ind))
   res$X.use <- X[XX.ind]
   res$Y.use <- Y[XX.ind]

   res$order.old<-order(c(XX.ind, X0.ind, X1.ind))
   res$n.samp <- length(res$Y.use)	 
   res$d <- cbind(res$X.use, res$Y.use)

   ## check survey data
   if (any(supplement <0) || any(supplement >1)) 
      stop("survey data have to be between 0 and 1.")
   if(is.null(supplement))
     res$survey.samp <- res$survey.data <- res$survey.yes <- 0
   else
     if (dim(supplement)[2] != ndim)
       stop("when context=TRUE, use n by 3. Otherwise use n by 2 matrix for survey data")
     else {
       res$survey.samp <- length(supplement[,1])
       res$survey.data <- as.matrix(supplement)
       res$survey.yes <- 1
     }

   return(res)

}
