# A. New function: binary logit + maximum likelihood + optim()
fitBina <- function(data, name.y, name.x, par0 = NULL, 
  model = c("logistic", "probit"), ...) {
  # A1. Inputs
  var.y <- data[, name.y, drop = FALSE]
  var.x <- cbind(Intercept = 1, data[, name.x])
  y <- as.matrix(var.y); x <- as.matrix(var.x)
  if (is.null(par0)) {par0 <- rep(x = 0, times = ncol(x))}
  model <- match.arg(model)
  cdf <- switch(EXPR = model, logistic = plogis, probit = pnorm)

  # A2. Construct loglikehood function and maximization by optim()
  loglike <- function(b, x, y) {
    p <- as.vector(cdf(x %*% b))
    val <- sum(y * log(p) + (1 - y) * log(1 - p))
    return(val)
  }
  best <- optim(par = par0, fn = loglike, x = x, y = y, hessian = TRUE,
    control = list(fnscale = -1), ...)

  # A3. Compute coeffient and p value
  coeff <- best$par
  covar <- solve(best$hessian)
  error <- matrix(data = diag(-1 * covar) ^ 0.5, ncol = 1)
  zvalu <- coeff / error
  pvalu <- 2 * pt(abs(zvalu), df = nrow(x) - ncol(x), lower.tail = FALSE)
  out <- data.frame(round(cbind(coeff, error, zvalu, pvalu), 3))
  colnames(out) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(out) <- colnames(x)
  out$sig <- ifelse(test = out[, 4] >= 0.10, yes = ' ', 
        no = ifelse(test = out[, 4] >= 0.05, yes = '.',    
        no = ifelse(test = out[, 4] >= 0.01, yes = '*', 
        no = ifelse(test = out[, 4] >= 0.001, yes = '**', no = '***'))))

  # A4. Output
  result <- list(data = data, name.y = name.y, name.x = name.x,
    best = best, out = out, call = sys.call())
  class(result) <- "fitBina"
  return(result)
}
print.fitBina <- function(x, ...) {print(x$out)}  # print method

# B. Test on real data with different initial value, methods, and models
library(erer); data(daIns)
nx <- colnames(daIns)[-c(1, 14)]
ha <- fitBina(data = daIns, name.y = "Y", name.x = nx, method = "BFGS")
hb <- update(ha, model = "probit")
hc <- update(ha, par0 = c(-4, 0.2, 0.01, 0.8, 0.1, -0.3, -0.3, 0.01,
  1.6, -0.2, -0.01, 0, 0))
hd <- update(ha, method = "Nelder-Mead")
  
names(ha); names(ha$best)
rbind(ha$best$counts, hb$best$counts, hc$best$counts, hd$best$counts)
ha$out[1:3, ]  # my binary logit model
hb$out[1:3, ]  # my binary probit model
                                                    
# C. Comparison with base R
ka <- glm(formula = Y ~ Injury + HuntYrs + Nonres + Lspman + Lnong +
  Gender + Age + Race + Marital + Edu + Inc + TownPop,
  family = binomial(link = "logit"), data = daIns, x = TRUE)
kb <- update(ka, family = binomial(link = "probit"))
round(x = summary(ka)$coefficients, digits =3)[1:3, ]  # R logit
round(x = summary(kb)$coefficients, digits =3)[1:3, ]  # R pro