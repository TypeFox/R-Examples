meanp <-
function(p) {
   keep <- (p >= 0) & (p <= 1)
   pi <- mean(p[keep])
   k <- length(p[keep])
   if(sum(1L * keep) < 4)
      stop("Must have at least four valid p values")
   z <- (0.5 - pi) * sqrt(12 * k)
   if(k != length(p)) {
      warning("Some studies omitted")
   }
   res <- list(z = z, p = pnorm(z, lower.tail = FALSE),
      validp = p[keep])
   class(res) <- c("meanp", "metap")
   res
}
print.meanp <- function(x, ...) {
   cat("z = ", x$z, " p = ", x$p, "\n")
   invisible(x)
}
