`ojaRank` <- function(X, x = NULL, p = NULL, silent = FALSE, na.action = na.fail){
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

   p <- checkP(p = p,n = n, k = k, silent = silent, type = "rank")
   
   if (!is.null(x)){
      if (p < 1) subSampleMessage(p, silent, string = "rank")
      rnk <- ojaGradient.hyperplanes(RCM.hyperplanes(X, p = p), x)
      return(rnk)
   }else{
       if ((p >= 1) && (k <= 10)){
         rk <- rep(1, n*k)
         out <- .C("ojacrnk", as.double(c(t(X))), as.integer(n), as.integer(k), ans = as.double(rk))
         OjaRanks <- matrix(out$ans, n, k, byrow=T)
      }else{
         subSampleMessage(p, silent, string = "ranks")
         Hyperplanes <- RCM.hyperplanes(X, p = p)      
         OjaRanks <- apply(X, 1, function(y){ojaGradient.hyperplanes(Hyperplanes, y)})
         if (is.matrix(OjaRanks)){
            OjaRanks <- t(OjaRanks)
         }else{
            OjaRanks <- as.matrix(OjaRanks)
         }
      }
      rownames(OjaRanks) <- rownames(X)
      colnames(OjaRanks) <- NULL
      return(OjaRanks)
   }
}
