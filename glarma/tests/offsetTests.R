### model.frame
### glm example
require(MASS)
data(Insurance)
glmmod <- glm(Claims ~ District + Group + Age + offset(log(Holders)),
              data = Insurance, family = poisson)
head(model.frame(glmmod))
### glarma example
require(glarma)
data(DriverDeaths)
y <- DriverDeaths[, "Deaths"]
X <- as.matrix(DriverDeaths[, 2:5])
Population <- DriverDeaths[, "Population"]

### No ARMA component
glarmamodNoARMA <- glarma(y, X, offset = log(Population/100000),
                          type = "Poi", method = "FS",
                          residuals = "Pearson", maxit = 100, grad = 1e-6)
head(model.frame(glarmamodNoARMA))
glmmod <- glm(y ~ X - 1, offset = log(Population/100000),
              family = poisson)
head(model.frame(glmmod))
summary(glarmamodNoARMA)
summary(glmmod)
print(glarmamodNoARMA)
print(glmmod)

### No offset
glarmamod <- glarma(y, X, phiLags = c(12),
                    type = "Poi", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)
head(model.frame(glarmamod))

### Offset included
glarmamodOffset <- glarma(y, X, offset = log(Population/100000),
                          phiLags = c(12),
                          type = "Poi", method = "FS",
                          residuals = "Pearson", maxit = 100, grad = 1e-6)
head(model.frame(glarmamodOffset))

### summary
summary(glmmod)
summary(glarmamodOffset)
summary(glarmamod)

### print
print(glmmod)
print(glarmamodOffset)
print(glarmamod)

### coef
coef(glmmod)
coef(glarmamodOffset)
coef(glarmamod)





