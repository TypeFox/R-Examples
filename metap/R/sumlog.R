sumlog <-
function(p) {
   keep <- (p > 0) & (p <= 1)
   lnp <- log(p[keep])
   chisq <- (-2) * sum(lnp)
   df <- 2 * length(lnp)
   if(sum(1L * keep) < 2)
      stop("Must have at least two valid p values")
   if(length(lnp) != length(p)) {
      warning("Some studies omitted")
   }
   res <- list(chisq = chisq, df = df,
      p = pchisq(chisq, df, lower.tail = FALSE), validp = p[keep])
   class(res) <- c("sumlog", "metap")
   res
}
print.sumlog <- function(x, ...) {
   cat("chisq = ", x$chisq, " with df = ", x$df, " p = ", x$p, "\n")
   invisible(x)
}
