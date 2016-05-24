predict.lqa <-
function (object, new.x = NULL, new.y = NULL, weights = rep (1, n.newobs), ...)
{
  if (is.null (new.x))
    new.x <- object$x

  new.x <- as.matrix (new.x)

  if (dim (new.x)[2] == 1)   # in this case new.x was a vector and hence we must transpose it
    new.x <- t (new.x)

  if (dim (new.x)[2] != length (object$coef))
    stop ("'new.x' does not match number of regressors in 'object'")

  n.newobs <- nrow (new.x)  
  new.dev <- NULL

  family <- object$family
  coef <- object$coefficients
  eta.new <- drop (new.x %*% coef)
  mu.new <- family$linkinv (eta.new)

  rownames.x <- if (is.null (rownames (new.x)))
                   1 : dim (new.x)[1]
                else
                   rownames (new.x)

  names (eta.new) <- names (mu.new) <- rownames.x

  if (!is.null (new.y))
    new.dev <- sum (family$dev.resids (new.y, mu.new, weights))

  pred.obj <- list (deviance = new.dev, tr.H = object$rank, n.newobs = n.newobs, eta.new = eta.new, mu.new = mu.new, lqa.obj = object, new.y = new.y)
  class (pred.obj) <- "pred.lqa"  
  pred.obj
}

