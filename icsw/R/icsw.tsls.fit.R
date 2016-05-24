icsw.tsls.fit <- function(
  D, X, Y, Z, W,
  weights, estimand = c("ATE", "ATT"),
  min.prob.quantile = NULL,
  min.prob = NULL, ...
  ) {

  complier.fit <- compliance.score(D, Z, W, weights = weights, ...)
  scores <- if (estimand[1] == "ATE") {
    complier.fit$C.score
  } else if (estimand[1] == "ATT") {
    (complier.fit$C.score /
     (complier.fit$A.score + complier.fit$C.score))
  } else {
    stop("Invalid estimand.")
  }

  if (!is.null(min.prob.quantile))
    min.prob <- quantile(scores[scores > 0], min.prob.quantile)

  tsls.wfit(
    cbind(X, D), Y, cbind(X, Z),
    weights = weights / clip.small.probs(scores, min.prob)
    )
}
