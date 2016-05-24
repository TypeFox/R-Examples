votep <-
function(p, alpha = c(0.5, 0.5)) {
   alpha <- ifelse(alpha > 1, alpha / 100, alpha) # if percent
   stopifnot(alpha > 0, alpha < 1)
   stopifnot(alpha[1] <= alpha[2])
   keep <- (p >= 0) & (p <= 1)
   pi <- p[keep]
   k <- length(pi)
   if(sum(1L * keep) < 2)
      stop("Must have at least two valid p values")
   pos <- sum(1L * (pi < alpha[1]))
   neg <- sum(1L * (pi > alpha[2]))
   if(k != length(p)) {
      warning("Some studies omitted")
   }
   if((pos + neg) <= 0) {
      warning("All p values are within specified limits of alpha")
      p <- 1
   } else {
      p = binom.test(pos, pos + neg, 0.5, alternative = "greater")$p.value
   }
   res = list(p = p, pos = pos, neg = neg, alpha = alpha, validp = pi)
   class(res) <- c("votep", "metap")
   res
}
print.votep <- function(x, ...) {
   cat("p = ", x$p, "\n")
   invisible(x)
}
