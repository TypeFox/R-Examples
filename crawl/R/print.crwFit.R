#' @method print crwFit
#' @export
print.crwFit <- function(x, ...)
{
  fit <- x
  cat("\n\n")
  cat("Continuous-Time Correlated Random Walk fit\n\n")
  cat('Models:\n')
  cat("--------\n")
  cat("Movement   "); cat(as.character(fit$mov.model)); cat("\n")
  cat("Error   "); cat(as.character(fit$err.model)); cat("\n")
  if (fit$random.drift | !is.null(fit$stop.model)) cat("with ")
  if (fit$random.drift) cat("Random Drift")
  if (fit$random.drift & !is.null(fit$stop.model)) cat(" and ")
  if (!is.null(fit$stop.model)) cat("Movement Stops")
  cat("\n\n")
  out <- as.data.frame(round(cbind(fit$par, fit$se, fit$ci[, 1], fit$ci[, 2]), 3))
  colnames(out) <- c("Parameter Est.", "St. Err.", "95% Lower", "95% Upper")
  rownames(out) <- fit$nms
  out[!is.na(fit$fixPar), 2:4] <- "."
  print(out)
  cat("\n\n")
  cat(paste("Log Likelihood =", round(fit$loglik, 3),"\n", collapse=""))
  cat(paste("AIC =", round(fit$aic, 3),"\n", collapse=""))
  cat("\n\n\n")
}

