SSallometric <-
structure(function (x, alpha, beta) 
{
    .expr1 <- x^beta
    .expr3 <- log(x)
    .expr4 <- .expr1 * .expr3
    .value <- alpha * .expr1
    .grad <- array(0, c(length(.value), 2L), list(NULL, c("alpha", 
        "beta")))
    .hessian <- array(0, c(length(.value), 2L, 2L), list(NULL, 
        c("alpha", "beta"), c("alpha", "beta")))
    .grad[, "alpha"] <- .expr1
    .hessian[, "alpha", "alpha"] <- 0
    .hessian[, "alpha", "beta"] <- .hessian[, "beta", "alpha"] <- .expr4
    .grad[, "beta"] <- alpha * .expr4
    .hessian[, "beta", "beta"] <- alpha * (.expr4 * .expr3)
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
}, initial = function (mCall, data, LHS) {
  xy <- data.frame(sortedXyData(mCall[["x"]], LHS, data))
  if (nrow(xy) < 3) 
    stop("Too few observations to fit allometric function")
  pars <- as.vector(coef(lm(I(log(y)) ~ I(log(x)), 
                            data = xy)))
  pars[1] <- exp(pars[1])
  names(pars) <- mCall[c("alpha", "beta")]
  return(pars)
}
, pnames = c("alpha", "beta"), class = "selfStart")
