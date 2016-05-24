asGnm.unidiff <- asGnm.assocmod <- function(object, ...) {
    object$call <- object$call.gnm
    class(object) <- c("gnm", "glm", "lm")
    object
}

printModelHeading <- function(x, digits=max(3, getOption("digits") - 4)) {
  cat("Call:\n", deparse(x$call), "\n", sep="", fill=TRUE)

  cat("Deviance Residuals:\n")
  if (x$df.residual > 5) {
      x$deviance.resid <- quantile(x$deviance.resid, na.rm = TRUE)
      names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q", "Max")
  }
  print(x$deviance.resid, digits=digits, na.print="", print.gap=2)
}

printModelStats <- function(x, digits=max(3, getOption("digits") - 4)) {
  cat("\nDeviance:              ", format(x$deviance, digits),
      "\nPearson chi-squared:   ",
      format(sum(na.omit(c(residuals(x, type="pearson")))^2), digits),
      "\nDissimilarity index:   ",
      format(sum(na.omit(c(abs(residuals(x, "response")))))/sum(na.omit(c(abs(fitted(x)))))/2*100, digits), "%",
      "\nResidual df:           ", x$df.residual,
      "\nBIC:                   ", x$deviance - log(sum(na.omit(c(x$data)))) * x$df.residual,
      "\nAIC:                   ", x$deviance - 2 * x$df.residual, "\n", sep="")
}

get.probs.asymm <- function(weighting, weights, indices) {
  # Weight with marginal frequencies, cf. Clogg & Shihadeh (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  if(weighting == "marginal")
      prop.table(rowSums(weights[indices,, drop=FALSE]))
  else if(weighting == "uniform")
      rep(1/nrow(weights), nrow(weights))
  else
      rep(1, nrow(weights))
}

get.probs.symm <- function(weighting, row.weights, col.weights, indices) {
  # Weight with marginal frequencies, cf. Clogg & Shihadeh (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  if(weighting == "marginal")
      prop.table(rowSums(row.weights[indices,, drop=FALSE]) + rowSums(col.weights[indices,, drop=FALSE]))
  else if(weighting == "uniform")
      rep(1/nrow(row.weights), nrow(row.weights))
  else
      rep(1, nrow(row.weights))
}

get.probs <- function(ass) {
  if(inherits(ass, "assoc.symm")) {
      if(length(ass$row.sup) > 0 & length(ass$col.sup) > 0) {
          p <- c(get.probs.symm(ass$weighting, ass$row.weights, ass$col.weights, setdiff(seq(nrow(ass$row)), ass$row.sup)),
                 get.probs.symm(ass$weighting, ass$row.weights, ass$col.weights, ass$row.sup))
          list(rp=p, cp=p)
      }
      else {
          p <- get.probs.symm(ass$weighting, ass$row.weights, ass$col.weights, seq(nrow(ass$row)))
          list(rp=p, cp=p)
      }

  }
  else {
      if(length(ass$row.sup) > 0)
          rp <- c(get.probs.asymm(ass$weighting, ass$row.weights, setdiff(seq(nrow(ass$row)), ass$row.sup)),
                  get.probs.asymm(ass$weighting, ass$row.weights, ass$row.sup))
      else
          rp <- get.probs.asymm(ass$weighting, ass$row.weights)

      if(length(ass$col.sup) > 0)
          cp <- c(get.probs.asymm(ass$weighting, ass$col.weights[setdiff(seq(ncol(ass$col)), ass$col.sup),]),
                  get.probs.asymm(ass$weighting, ass$col.weights[ass$col.sup,]))
      else
          cp <- get.probs.asymm(ass$weighting, ass$col.weights)

      list(rp=rp, cp=cp)
  }
}
