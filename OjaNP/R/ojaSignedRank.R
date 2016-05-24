`ojaSignedRank` <- function(X, x = NULL, p = NULL, silent = FALSE, na.action = na.fail){
   X <- na.action(X)
   if (!all(sapply(X, is.numeric))) stop("'X' must be numeric")
   X <- as.matrix(X)
   k <- ncol(X)
   n <- nrow(X)

   if (!is.null(x)){
      if (!is.numeric(x)) stop("'x' must be numeric")
      if (length(x) != k) stop("'x' and 'X' must have the same dimension")
      x <- as.vector(x)
   }
   
   p <- checkP(p = p,n = n, k = k, silent = silent, type = "signedrank")

   if (!is.null(x)){
      if (p < 1) subSampleMessage(p, silent, string = "signed rank")
      OjaSignedRanks <- ojaGradient.hyperplanes(RCM.hyperplanes3(X, p = p), x)
      return(OjaSignedRanks)
   }else{
       if (p >= 1){
         rk <- rep(1, n*k)
         out <- .C("ojasrnk", as.double(c(t(X))), as.integer(n), as.integer(k), ans = as.double(rk))
         OjaSignedRanks <- matrix(out$ans, n, k, byrow=T)
      }else{
         subSampleMessage(p, silent, string = "signed ranks")
         Hyperplanes <- RCM.hyperplanes3(X, p = p)
         OjaSignedRanks <- apply(X, 1, function(y){ojaGradient.hyperplanes(Hyperplanes, y)})
         if (is.matrix(OjaSignedRanks)){
            OjaSignedRanks <- t(OjaSignedRanks)
         }else{
            OjaSignedRanks <- as.matrix(OjaSignedRanks)
         }
      }
      rownames(OjaSignedRanks) <- rownames(X)
      colnames(OjaSignedRanks) <- NULL
      return(OjaSignedRanks)
   }
}
