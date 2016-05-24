`ojaSign` <- function(X, x = NULL, center = "ojaMedian", p = NULL, silent = FALSE, na.action = na.fail, ...){
   X <- na.action(X)
   if (!all(sapply(X, is.numeric)))  stop("'X' must be numeric")
   X <- as.matrix(X)   
   k <- ncol(X)
   n <- nrow(X)
   if (k<=1) stop("Dimension must be at least 2. For univariate signs use 'sign()'.")
   
   if (!is.null(x)){
      if (!is.numeric(x)) stop("'x' must be numeric")
      if (length(x) != k) stop("'x' and 'X' must have the same dimension")
      x <- as.vector(x)
   } 
       
   CENTER <- checkCenter(center = center, X = X,...)
   p <- checkP(p = p, n = n, k = k, silent = silent, type = "sign")
   
   if (!is.null(x)){
      if (p < 1) subSampleMessage(p,silent,string = "sign")
      sgn <- ojaGradient.hyperplanes(SCM.hyperplanes(X, center = CENTER, p = p, ...), x)
      return(sgn)
   }else{
      if ((p >= 1) && (k <= 10)){
         X <- sweep(X, 2, CENTER)
         rk <- rep(1, n*k)
         out <- .C("ojasn", as.double(c(t(X))), as.integer(n), as.integer(k), ans = as.double(rk))
         OjaSigns <- matrix(out$ans, n, k, byrow=T)
      }else{
         subSampleMessage(p, silent, string = "signs")
         Hyperplanes <- SCM.hyperplanes(X, center = CENTER, p = p, ...)      
         OjaSigns <- t(apply(X, 1, function(y){ojaGradient.hyperplanes(Hyperplanes, y)}))
      }
      rownames(OjaSigns) <- rownames(X)
      colnames(OjaSigns) <- NULL
      return(OjaSigns)
   }
}
