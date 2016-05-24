#
# rv package
# Simple Regression imputation example
# 2007-07-05
#

require(rv)

setnsims(1000)

y <- c(67, 100, 45, 69, 80, 74, 76, 89, 72, 66, NA, NA, NA, NA, NA)
x <- c(85,  96, 75, 56, 45, 75, 95, 96, 98, 80, 54, 79, 42, 67, 86)
y.missing <- is.na(y)

missing.students <- c('Alice','Bob','Cecil','Dave','Ellen')

lm.0 <- lm(y ~ x) # na.omit!

par(mfrow=c(2,2))

plot(x, y, xlim=c(0,100),ylim=c(0,100))
title('Data and fitted regression line')
abline(lm.0, col='gray')


post <- postsim(lm.0)

sigma <- post['sigma']
beta  <- post[-1]

cat("Posterior simulations of the linear regression fit:\n")
print(sigma)
print(beta)

X.missing <- cbind(1, x[is.na(y)])
y.pred <- rvnorm( mean=X.missing %*% beta, sd=sigma )

names(y.pred) <- missing.students
y[is.na(y)] <- y.pred

cat("Predicted scores:\n")
print(y.pred)

cat("Imputed scores:\n")
print(y)

cat("Class median:\n")
print(median.rv(y))

cat("Class mean:\n")
print(mean(y))


plot(x, y,xlim=c(0,100),ylim=c(0,100),rvcol="red")
title('20 posterior simulations of the regression line')
invisible(rvpar(line.sample=20))
abline(beta, col='gray')

plot(x, y, xlim=c(0,100),ylim=c(0,100), xlab="midterm", ylab="final", rvcol="red")
title(main="Intervals for predicted examination scores")
abline(lm.0, col='gray')

y.med <- rvmedian(y)
plot(x, y.med, xlim=c(0,100),ylim=c(0,100), ylab="y (imputed: solid circles)")
title(main="Imputed posterior medians")
abline(lm.0, col='gray')
points(x[y.missing], y.med[y.missing], pch=18)
