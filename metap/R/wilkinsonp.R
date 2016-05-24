wilkinsonp <-
function(p, r = 1, alpha = 0.05) {
   alpha <- ifelse(alpha > 1, alpha / 100, alpha) # if percent
   stopifnot(alpha > 0, alpha < 1)
   alpha <- ifelse(alpha > 0.5, 1 - alpha, alpha) # map to 0 to 0.5
   keep <- (p >= 0) & (p <= 1)
   pi <- p[keep]
   k <- length(pi)
   if(sum(1L * keep) < 2)
      stop("Must have at least two valid p values")
   if(k != length(p)) {
      warning("Some studies omitted")
   }
   if((r < 1) | (r > k)) {
      r <- 1
      warning("Illegal r set to 1")
   }
   pi <- sort(pi)
   pr <- pi[r]
   res <- list(p = pbeta(pr, r, k + 1 - r), pr = pr,
      r = r, critp = qbeta(alpha, r, k + 1 - r), alpha = alpha, validp = pi)
   class(res) <- c("wilkinsonp", "metap")
   res
}
minimump <- function(p,  alpha = 0.05) {
   res <- wilkinsonp(p, r = 1, alpha)
   class(res) <- c("minimump", "metap")
   res
}
print.wilkinsonp <- function(x, ...) {
   r <- x$r
   p <- x$p
   cat("p = ", p, " using ", paste(r, "th minimum p", sep = ""), "\n")
   invisible(x)
}
print.minimump <- function(x, ...) {
   r <- x$r
   p <- x$p
   cat("p = ", p, " using minimum p", "\n")
   invisible(x)
}
