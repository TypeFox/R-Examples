## Example of prediction variance computations

## Predictions
x <- rnorm(15)
x <- sort(x)
y <- x + rnorm(15)
predict(lm(y ~ x))
new <- data.frame(x = x)
predict(lm(y ~ x), new, se.fit = TRUE)
pred.w.plim <- predict(lm(y ~ x), new, interval="prediction")
pred.w.clim <- predict(lm(y ~ x), new, interval="confidence")

model <- regress(y~x,verbose=10)
##plot(x,y)
lm(y~x)$df.residual
two <- qt((1-0.95)/2,13)
hwid <-  with(model,two*c(1,-1) %o% sqrt(predictedVariance2))
pred.w.clim2 <- with(model,cbind(fitted,t(hwid)+as.vector(fitted)))
hwid <-  with(model,two*c(1,-1) %o% sqrt(predictedVariance2+predictedVariance))
pred.w.plim2 <- with(model,cbind(fitted,t(hwid)+as.vector(fitted)))

## Test if predictions and prediction intervals are the same
if(sum((pred.w.plim - pred.w.plim2)^2) + sum((pred.w.clim - pred.w.clim2^2)) > 1e-10) stop("Error with prediction variances")

par(mfrow=c(1,1))
matplot(new$x,cbind(pred.w.clim, pred.w.plim[,-1]),
        lty=c(1,2,2,3,3), type="l", ylab="predicted y",lwd=3)
points(x,y)

lines(x,model$fitted)
##lines(x,model$fitted-two*sqrt(model$predictedVariance))
##lines(x,model$fitted+two*sqrt(model$predictedVariance))
lines(x,model$fitted-two*sqrt(model$predictedVariance2))
lines(x,model$fitted+two*sqrt(model$predictedVariance2))
lines(x,model$fitted-two*sqrt(model$predictedVariance+model$predictedVariance2))
lines(x,model$fitted+two*sqrt(model$predictedVariance+model$predictedVariance2))

## Test if predictions and prediction intervals are the same
if(sum((pred.w.plim - pred.w.plim2)^2) + sum((pred.w.clim - pred.w.clim2^2)) > 1e-10) stop("Error with prediction variances")

