library(spls)
data(yeast)

set.seed(1)
cv <- cv.spls(yeast$x, yeast$y, eta = seq(0.1,0.9,0.1), K = c(5:10))
print(cv)

## Perform SPLS with eta=0.7, 8 latent components
f <- spls(yeast$x, yeast$y, K=8, eta=0.7)
print(f)

## Extract and print out coefficients
coef.f <- coef(f)
coef.f[, 1]

predict.f <- predict(f, type="fit")
predict.f

## Coefficient path plot
plot(f, yvar=1)
dev.new()

## Coefficient plot of selected variables
coefplot.spls(f, xvar=c(1:4))
