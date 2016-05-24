# A. Define a class
setClass(Class = "lm2", slots = c(data = "data.frame", name.y ="character",
  name.x = "character", fitted = "matrix", resid = "matrix", 
  sigma = "numeric", covar = "matrix", out = "ANY")
)

# B. Object constructor 
lm2 <- function(data, name.y, name.x) { 
  # B1. Examine inputs
  if (!inherits(data, "data.frame")) stop("y should be data.frame.\n")
  if (!inherits(name.y, "character") || !inherits(name.x, "character")) {
    stop("name.y and name.x should be character.\n")
  }
  if (!all(c(name.y, name.x) %in% colnames(data))) {
    stop("name.y or name.x is not in the data.\n")
  }    
  var.y <- data[, name.y, drop = FALSE]
  var.x <- cbind(Intercept = 1, data[, name.x])
  y <- as.matrix(var.y); x <- as.matrix(var.x)
  
  # B2. Compute coeff, fitted value, residuals, sigma, var, t, p
  coeff <- solve(t(x) %*% x) %*% t(x) %*% y
  fitted <- x %*% coeff
  resid <- y - fitted
  sigma2 <- crossprod(resid) / (nrow(x) - nrow(coeff))
  sigma <- c(sqrt(sigma2))
  covar <- c(sigma2) * solve(crossprod(x))
  error <- matrix(data = diag(covar)^0.5, ncol = 1)
  tvalu <- coeff / error
  pvalu <- 2 * pt(abs(tvalu), df = nrow(x) - ncol(x), lower.tail=FALSE)
  
  # B3. Outputs
  out <- cbind(round(x = coeff, digits = 7), round(x = error, digits = 7), 
               round(x = tvalu, digits = 3), round(x = pvalu, digits = 6))
  colnames(out) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  result <- new(Class = "lm2", data = data, name.y = name.y, 
    name.x = name.x, fitted = fitted, resid = resid, sigma = sigma, 
    covar = covar, out = out)
  return(result)
}

# C. Define methods 
setMethod(f = "show", signature = signature(object = "lm2"), 
  definition = function(object) {
    cat("\n________________ My OLS Results ________________\n")
    print(object@out)
    invisible(object)
  }
)

setMethod(f = "plot", signature = signature(x = "lm2"),
  definition = function(x) {
    plot(x@resid, main = "Residuals from my OLS regression")
    invisible(NULL)
  }
)