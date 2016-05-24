 structSumIPC <- function(x, ord = c("PA", "BC", "DE", "FG", "HI", "JK", "LM", "NO")) {
  if((length(x) / 8) != round(length(x) / 8)) {stop("length of x must be a mutliple of 8.")}
  if(!is.data.frame(x)) {
    mat <- matrix(x, ncol=8, byrow=FALSE)
    x <- data.frame(mat)
  }
  colnames(x) <- ord

  DOM <- .25 * (x$PA - x$HI + .71 * (x$NO + x$BC - x$FG - x$JK))
  LOV <- .25 * (x$LM - x$DE + .71 * (x$NO - x$BC - x$FG + x$JK))
  DEG <- atan2(DOM, LOV) * (180 / pi)
  DEG <- ifelse(DEG >= 0, DEG, DEG + 360)
  AMP <- sqrt(DOM^2 + LOV^2)
  ELEV <- rowMeans(x)
  SStot <- rowSums((x - matrix(ELEV, nrow=nrow(x), ncol=8, byrow=FALSE))^2)
  Rsq <- 4 * AMP^2 / SStot

  res <- data.frame(DOM, LOV, DEG, AMP, ELEV, SStot, Rsq)
  res
}
