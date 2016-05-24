Kappa <- function (x, weights = c("Equal-Spacing", "Fleiss-Cohen"))
{
  if (is.character(weights))
      weights <- match.arg(weights)

  d  <- diag(x)
  n  <- sum(x)
  nc <- ncol(x)
  colFreqs <- colSums(x)/n
  rowFreqs <- rowSums(x)/n

  ## Kappa
  kappa <- function (po, pc)
    (po - pc) / (1 - pc)
  std  <- function (p, pc, kw, W = diag(1, ncol = nc, nrow = nc)) {
    sqrt((sum(p * sweep(sweep(W, 1, W %*% colSums(p) * (1 - kw)), 2, W %*% rowSums(p) * (1 - kw)) ^ 2) - (kw - pc * (1 - kw)) ^ 2) / crossprod(1 - pc) / n)
  }

  ## unweighted
  po <- sum(d) / n
  pc <- crossprod(colFreqs, rowFreqs)[1]
  k <- kappa(po, pc)
  s <- std(x / n, pc, k)

  ## weighted
  W <- if (is.matrix(weights))
    weights
  else if (weights == "Equal-Spacing")
    1 - abs(outer(1:nc, 1:nc, "-")) / (nc - 1)
  else
    1 - (abs(outer(1:nc, 1:nc, "-")) / (nc - 1))^2
  pow <- sum(W * x) / n
  pcw <- sum(W * colFreqs %o% rowFreqs)
  kw <- kappa(pow, pcw)
  sw <- std(x / n, pcw, kw, W)

  structure(
            list(Unweighted = c(
                   value = k,
                   ASE   = s
                   ),
                 Weighted = c(
                   value = kw,
                   ASE   = sw
                   ),
                 Weights = W
                 ),
            class = "Kappa"
       )
}

print.Kappa <-
		function (x, digits=max(getOption("digits") - 3, 3), CI=FALSE, level=0.95, ...)
{
	tab <- rbind(x$Unweighted, x$Weighted)
	z <- tab[,1] / tab[,2]
	tab <- cbind(tab, z, `Pr(>|z|)` = 2 * pnorm(-abs(z)))
	if (CI) {
		q <-  qnorm((1 + level)/2)
		lower <- tab[,1] - q * tab[,2]
		upper <- tab[,1] + q * tab[,2]
		tab <- cbind(tab, lower, upper)
	}

	rownames(tab) <- names(x)[1:2]
	print(tab, digits=digits, ...)
	invisible(x)
}

summary.Kappa <- function (object, ...)
  structure(object, class = "summary.Kappa")

print.summary.Kappa <- function (x, ...) {
  print.Kappa(x, ...)
  cat("\nWeights:\n")
  print(x$Weights, ...)
  invisible(x)
}

confint.Kappa <- function(object, parm, level = 0.95, ...) {
  q <- qnorm((1 + level) / 2)
  matrix(c(max(-1, object[[1]][1] - object[[1]][2] * q),
           min(1, object[[1]][1] + object[[1]][2] * q),
           max(-1, object[[2]][1] - object[[2]][2] * q),
           min(1, object[[2]][1] + object[[2]][2] * q)),
         ncol = 2, byrow = TRUE,
         dimnames = list(Kappa = c("Unweighted","Weighted"), c("lwr","upr"))
         )
}

