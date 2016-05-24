.max.bias <- function(resi, x, ng = 10) {
	t2 <- range(x)
	MB <- tapply(resi, cut(x, seq(from = t2[1], to = t2[2], length.out = ng+1), include.lowest=TRUE), mean)
  max(abs(MB), na.rm = TRUE)
}

.rmse <- function(x) {
  sqrt(mean(x^2, na.rm=TRUE))
}

.r2 <- function(x, obs) {
  cor(x, obs, use="pairwise.complete.obs")^2
}

.get.rand <- function(range=NULL) {
   x <- as.double(0)
   x <- .C("GetSetRand", as.vector(x), PACKAGE="rioja")[[1]]
   if (is.null(range))
      return (x)
   else
      return(as.integer(as.double(range) * x))
}

.set.rand.seed <- function(x=1) {
   if (x <= 0) 
     stop("Must set seed with x < 0")
   x <- .C("GetSetRand", as.double(-x), PACKAGE="rioja")
   for (i in 1:10)
     x <- .get.rand()
}

