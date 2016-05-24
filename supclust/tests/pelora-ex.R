library(supclust)
data(leukemia, package="supclust")

set.seed(724)
xN <- matrix(rnorm(750), 3, 250)
x  <- leukemia.x
dimnames(x) <- list(1:nrow(x), paste("V",1:ncol(x),sep=""))

.proctime00 <- proc.time()

## 1. Without clinical variables:
fit1 <- pelora(leukemia.x, leukemia.y, noc=3, lam=1/32, fl="pm", sta=FALSE)

## 2. With clinical variables:
fitW <- pelora(leukemia.x[,101:250],leukemia.y, leukemia.x[,1:100], noc = 3,
               flip = "pm", standardize = FALSE)

## 3. With dimnames 
fit3 <- pelora(x, leukemia.y, noc = 3, lambda = 1/32, flip="pm", stand = FALSE)

## 4. Test the tracing
fit4 <- pelora(x, leukemia.y, noc = 1, trace = 2, flip = "pm", stand = FALSE)

## 5. Without tracing
fit5 <- pelora(x, leukemia.y, noc = 10, trace = 0, flip = "cor", stand = FALSE)

## Running time
cat('Time elapsed: ', proc.time() - .proctime00,'\n')

## Checking the output of fit1
fit1
summary(fit1)
coef(fit1)
plot(fit1)
fitted(fit1)

identical(predict(fit1), fitted(fit1))
predict(fit1, type = "cla")
predict(fit1, type = "pro")

predict(fit1, newdata = xN)
predict(fit1, newdata = xN, type = "pro")
predict(fit1, newdata = xN, type = "cla")


## Checking the output of fit2
fitW
summary(fitW)
coef(fitW)
plot(fitW)
fitted(fitW)

identical(predict(fitW), fitted(fitW))
predict(fitW, type = "cla")
predict(fitW, type = "pro")

predict(fitW, newdata = xN[,101:250], newc = xN[,1:100])
predict(fitW, newdata = xN[,101:250], newc = xN[,1:100], ty = "pro")
predict(fitW, newdata = xN[,101:250], newc = xN[,1:100], ty = "cla", noc=c(1,3))

## Checking the output of fit3
fit3
summary(fit3)
coef(fit3)
plot(fit3)
fitted(fit3)

identical(predict(fit3), fitted(fit3))
predict(fit3, type = "cla")
predict(fit3, type = "pro")

predict(fit3, newdata = xN)
predict(fit3, newdata = xN, type = "pro")
predict(fit3, newdata = xN, type = "cla")

## checking the output of fit4
stopifnot(identical(fit1$genes[[1]], fit4$genes[[1]]))
str(fit4)
summary(fit4)

## Checking the output of fit3
fit5
summary(fit5)
coef(fit5)
plot(fit5)
fitted(fit5)

identical(predict(fit5), fitted(fit5))
predict(fit5, type = "cla")
predict(fit5, type = "pro")

predict(fit5, newdata = xN)
predict(fit5, newdata = xN, type = "pro")
predict(fit5, newdata = xN, type = "cla")

## Running time
cat('Time elapsed: (total) ', proc.time() - .proctime00,'\n')
