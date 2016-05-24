# S3method for "singh"
singh <-
function(data, cov, inverted = FALSE) UseMethod("singh")

# ---------------------------------
singh.default <-
function(data, cov, inverted = FALSE)
{
   if (!inherits(data, c("data.frame", "matrix")))
      stop("data must be a data.frame or matrix!")
   stopifnot(is.matrix(cov))
   if (ncol(data) != ncol(cov))
      stop("incompatible dimensions!")
   data <- as.matrix(data)
   n <- nrow(data)
   p <- nrow(cov)
   if (!inverted) {
      invS <- solve(cov)
   } else {
      invS <- cov
   }
   delta <- matrix(NA, choose(n, 2), p)
   l1 <- rep(1:n, each = n)
   l2 <- rep(1:n, times = n)
   aux <- rbind(l1[l1 < l2], l2[l1 < l2])
   for(i in 1:nrow(delta)) {
      delta[i, ] <- diff(data[aux[1:2, i], ])
   }
   stat <- matrix(NA, choose(n, 2), p)
   colnames(stat) <- colnames(data)
   for(j in 1:p) {
      for(i in 1:nrow(delta)) {
         stat[i, j] <- crossprod(delta[i, j] * delta[i, ], invS[, j])
      }
   }
   # output
   out1 <- abs(apply(stat, 2, sum))
   ord <- order(out1, decreasing = TRUE)
   out1 <- out1[ord]
   out2 <- out1 / sum(out1)
   out3 <- cumsum(out2)
   out <- rbind(out1, out2, out3)
   rownames(out) <- c("Singh statistic", "Proportion",
      "Cumulative proportion")
   class(out) <- "singh"
   return(out)
}

# -----------------------------------
# plot method
plot.singh <- function(x, ...)
{
   lab <- paste(colnames(x), " (", 
      round(x[2, ], 3)*100, "%", ")", sep = "")
   pie(x[2, ], labels = lab,
      main = "Importance of variables", ...)
}
