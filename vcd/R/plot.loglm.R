plot.loglm <- function(x,
                       panel = mosaic,
                       type = c("observed", "expected"),
                       residuals_type = c("pearson", "deviance"),
		       gp = shading_hcl,
		       gp_args = list(),
                       ...)
{
  residuals_type <- match.arg(tolower(residuals_type), c("pearson", "deviance"))
  if(is.null(x$fitted)) x <- update(x, fitted = TRUE)
  expected <- fitted(x)
  residuals <- residuals(x, type = "pearson")
  observed <- residuals * sqrt(expected) + expected
  if(residuals_type == "deviance") residuals <- residuals(x, type = "deviance")

  gp <- if(inherits(gp, "grapcon_generator"))
    do.call("gp", c(list(observed, residuals, expected, x$df), as.list(gp_args))) else gp

  panel(observed, residuals = residuals, expected = expected, type = type,
    residuals_type = residuals_type, gp = gp, ...)
}

mosaic.loglm <- function(x, ...)
{
  plot(x, panel = mosaic, ...)
}

assoc.loglm <- function(x, ...)
{
  plot(x, panel = assoc, ...)
}

sieve.loglm <- function(x, ...)
{
  plot(x, panel = sieve, ...)
}
