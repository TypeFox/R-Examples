
init.vars <- function(model, debug = FALSE)
{
  if (is.null(model@sgfc)) {
    xreg <- stsm.sgf(model, gradient = FALSE, hessian = FALSE,
      deriv.transPars = FALSE)$constants
  } else
    xreg <- model@sgfc

  if (is.null(model@ssd)) {
    pg <- Mod(fft(model@diffy))^2 / (2 * pi * length(model@diffy))
  } else pg <- model@ssd

  y <- 2 * pi * pg

  fit <- lm.fit(x = xreg, y = y)
  vars <- coef(fit)

  if (debug)
  {
    stopifnot(all.equal(coef(lm(y ~ 0 + xreg)), vars, 
      check.attributes = FALSE))
  }

  if (!is.null(model@transPars))
  {
    msg <- paste("variance parameters are not reparameterized ", 
    "according to 'model@transPars' = '", model@transPars, "'.", sep = "")
    warning(msg)
  }

  return(list(vars = vars, fit = fit))
}
