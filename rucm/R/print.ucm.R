#'Print ucm Object
#'@method print ucm
#'@param x ucm object
#'@param ... Ignored.
#'@rdname print.ucm
#'@importFrom stats as.formula is.ts predict printCoefmat pt terms var
#'@export
print.ucm <- function(x, ...){
  var.est <- c(x$irr.var, x$est.var.level, x$est.var.slope, x$est.var.season, x$est.var.cycle)
  if (!is.null(x$est)){
  Estimate <- round(x$est,4)
  Approx.StdErr <- attr(x, "se")
  t.val <- Estimate/Approx.StdErr
  p.value = 2*pt(-abs(t.val), df = attr(x, "df"))
  TAB <- cbind(Estimate = Estimate, Approx.StdErr = Approx.StdErr, t.val = t.val, p.value = p.value)
  }
  cat("Call:\n")  
  print(x$call)
  cat("\nParameter estimates:\n")
  if (!is.null(x$est)){
    printCoefmat(TAB, P.values=TRUE, has.Pvalue=TRUE)
  } else print(x$est)
  cat("\nEstimated variance:\n")
  print(round(var.est,4), digits = 7)
}