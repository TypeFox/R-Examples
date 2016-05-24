sumz <-
function(p, weights = NULL, data = NULL, subset = NULL, na.action = na.fail)  {
   if(is.null(data)) data <- sys.frame(sys.parent())
   mf <- match.call()
   mf$data <- NULL
   mf$subset <- NULL
   mf$na.action <- NULL
   mf[[1]] <- as.name("data.frame")
   mf <- eval(mf, data)
   if(!is.null(subset)) mf <- mf[subset,]
   mf <- na.action(mf)
   p <- as.numeric(mf$p)
   weights <- mf$weights
   if(is.null(weights)) weights <- rep(1, length(p))
   if(length(p) != length(weights)) warning("Length of p and weights differ")
   keep <- (p > 0) & (p < 1)
   if(sum(1L * keep) < 2)
      stop("Must have at least two valid p values")
   if(sum(1L * keep) != length(p)) {
      warning("Some studies omitted")
      omitw <- weights[!keep]
      if(sum(1L * omitw) > 0) warning("Weights omitted too")
   }
   zp <- (qnorm(p[keep], lower.tail = FALSE) %*% weights[keep]) /
      sqrt(sum(weights[keep]^2))
   res <- list(z = zp, p = pnorm(zp, lower.tail = FALSE),
      validp = p[keep], weights = weights)
   class(res) <- c("sumz", "metap")
   res
}
print.sumz <- function(x, ...) {
   cat("sumz = ", x$z, "p = ", x$p, "\n")
   invisible(x)
}
