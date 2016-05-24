library("aroma.core")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Example 1: median polish
# From example(medpolish) in the 'stats' package
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Deaths from sport parachuting;  from ABC of EDA, p.224:
deaths <- matrix(c(14,15,14, 7,4,7, 8,2,10, 15,9,10, 0,2,0), ncol=3, byrow=TRUE)
rownames(deaths) <- c("1-24", "25-74", "75-199", "200++", "NA")
colnames(deaths) <- 1973:1975
print(deaths)

fit1 <- medpolish(deaths, trace=FALSE)
r1 <- residuals(fit1)
fit2 <- matrixBlockPolish(deaths, FUN=function(y, x, ...) median(y, ...))
r2 <- residuals(fit2)
stopifnot(all.equal(r1,r2))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Example 2: smooth spline polish ("spatial smoothing")
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# From example(image) in 'graphics' package
x <- y <- seq(-4*pi, 4*pi, len=27)
r <- sqrt(outer(x^2, y^2, FUN="+"))
z <- cos(r^2) * exp(-r/6)

fit <- matrixBlockPolish(z, FUN=function(z, x, ...) median(z, ...),
                                                      returnEffects=TRUE)
r1 <- residuals(fit)

fit <- matrixBlockPolish(z, FUN=function(z, x, ...) {
  fit <- smooth.spline(x=x, y=z, ...)
  predict(fit, x=x)$y
}, spar=0.25)
r2 <- residuals(fit)
print(range(r2))
image(r2)
