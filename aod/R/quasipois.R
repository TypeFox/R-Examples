quasipois <- function(formula, data, phi = NULL, tol = 0.001){
# check call validity
  CALL <- match.call(expand.dots = FALSE)
  mf <- model.frame(formula = formula, data = data)
  y <- model.response(mf)
# computations
  if(is.null(phi)){
    phi <- 1e-04
    X2 <- 0
    delta <- X2 + 2 * tol
    w <- rep(1, length(y))
    dfr <- data.frame(data, w = w)
    fm <- glm(formula = formula, family = poisson, weights = w, data = dfr)
    ok <- TRUE
    while(ok){
      X2 <- sum(residuals(fm, type = "pearson")^2)
      delta <- X2 - df.residual(fm)
      if(delta <= tol)
        ok <- FALSE
        else{
          phi <- phi * sum(residuals(fm, type = "pearson")^2) / df.residual(fm)
          w <- 1 / (1 + phi * fitted(fm))
          dfr <- data.frame(data, w = w)
          fm <- glm(formula = formula, family = poisson, weights = w, data = dfr)
          }
        }
      }
  else{
    fm <- glm(formula = formula, family = poisson, data = data)
    delta.logL <- 1
    while(delta.logL > 1e-06){
      w <- 1 / (1 + phi * fitted(fm))
      dfr <- data.frame(data, w = w)
      fm.new <- glm(formula = formula, family = poisson, weights = w, data = dfr)
      delta.logL <- logLik(fm.new) - logLik(fm)
      fm <- fm.new
      }
    }

# results
  w <- 1 / (1 + phi * fitted(fm))
  dfr <- data.frame(data, w = w)
  fm <- glm(formula = formula, family = poisson, weights = w, data = dfr)
# outputs
  new(Class = "glimQL", CALL = CALL, fm = fm, phi = phi)
  }

