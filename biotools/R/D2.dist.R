D2.dist <-
function(data, cov, inverted = FALSE)
{
   if (!inherits(data, c("data.frame", "matrix")))
      stop("data must be a data.frame or matrix!")
   stopifnot(is.matrix(cov))
   if (ncol(data) != ncol(cov))
      stop("incompatible dimensions!")
   x <- as.matrix(data)
   n <- nrow(x)
   D2 <- matrix(0, n, n)
   dimnames(D2) <- list(rownames(data), rownames(data))

   if (!inverted) {
      for(i in 1:n) {
         for(j in 1:n) {
         if (i > j) D2[i, j] <- crossprod((x[i, ] - x[j, ]),
    solve(cov, (x[i, ] - x[j, ])))
         }
      }
   } else {
     for(i in 1:n) {
         for(j in 1:n) {
         if (i > j) D2[i, j] <- crossprod((x[i, ] - x[j, ]),
    crossprod(cov, (x[i, ] - x[j, ])) )
         }
      }
   }
   return(as.dist(D2))
}
