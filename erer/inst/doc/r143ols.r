# A. Input data
library(erer); data(daIns)
in.data <- daIns
in.name.y <- "Y"
in.name.x <- colnames(daIns)[-c(1, 14)]; in.name.x

# B1. Raw data transformation
var.y <- in.data[, in.name.y, drop = FALSE]
var.x <- cbind(Intercept = 1, in.data[, in.name.x])
y <- as.matrix(var.y)
x <- as.matrix(var.x)

# B2. Compute coefficients, fitted value, residuals
coeff <- solve(t(x) %*% x) %*% t(x) %*% y
fitted <- x %*% coeff
resid <- y - fitted

# B3. Compute sigma, variance, t value, p value
sigma2 <- crossprod(resid) / (nrow(x) - nrow(coeff))
sigma <- sqrt(sigma2)
covar <- c(sigma2) * solve(crossprod(x))
error <- matrix(data = diag(covar) ^ 0.5, ncol = 1)
tvalu <- coeff / error
pvalu <- 2 * pt(abs(tvalu), df = nrow(x) - ncol(x), lower.tail = FALSE)

# C. Outputs
out <- cbind(round(x = coeff, digits = 7), round(x = error, digits = 7), 
             round(x = tvalu, digits = 3), round(x = pvalu, digits = 6))
colnames(out) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
sigma
out

# D. Comparison with 'lm'
gg <- lm(formula = Y ~ 1 + Injury + HuntYrs + Nonres + Lspman + Lnong +
  Gender + Age + Race + Marital + Edu + Inc + TownPop, data = daIns)
summary(gg)$r.squared; summary(gg)$sigma
summary(gg)