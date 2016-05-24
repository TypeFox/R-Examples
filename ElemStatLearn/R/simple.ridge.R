simple.ridge <- function(x, y, lambda=1, df, ...) {
   n <- NROW(x)
   p <- NCOL(x)
   mm <- colMeans(x)
   x <- scale(x, mm, TRUE)
   beta0 <- mean(y)  # intercept
   svd.x <- svd(x, nu=p, nv=p) # necessary?
   dd <- svd.x$d
   u  <- svd.x$u
   v  <- svd.x$v
   if (missing(df)) { df <- sapply(lambda, function(x) sum(dd^2/(dd^2+x)))
   } else {
   fun <- function(df,lambda) df-sum(dd^2/(dd^2+lambda))
   lambda <- sapply(df, FUN=function(df) uniroot(f=function(lambda) fun(df, lambda),
                     lower=-0.000001, upper=1000, maxiter=10000)$root)
   }# a good, general upper limit needs more thought
   # beta tiene que ser una matriz en el caso multiple, con
   # cada beta para cada lambda como columna:
   ystar <- t(u) %*% y
   dstar <- sapply(lambda, function(x) dd/(dd^2+x))
   beta <- matrix(NA, p, length(lambda))
   for (i in seq(along=lambda)) {
      beta[,i] <- sweep(v,2,dstar[,i], "*") %*% ystar }
   colnames(beta) <- lambda
   rownames(beta) <- colnames(x)
   return( list( beta0=beta0, beta=beta, lambda=lambda, df=df ) )
}

