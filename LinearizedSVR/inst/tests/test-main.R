
mse <- function(y, yhat) mean((y-yhat)^2)

## Check some basic results (this is really a classification problem)
set.seed(123)
test_that("Basic result should work", {
  dat <- rbind(data.frame(y=2, x1=rnorm(500, 1), x2=rnorm(500, 1)),
               data.frame(y=1, x1=rnorm(500,-1), x2=rnorm(500,-1)))
  mod <- LinearizedSVRTrain(X=as.matrix(dat[-1]), Y=dat$y, nump=6)
  res <- predict(mod, newdata=as.matrix(dat[-1]))
  expect_true(mean(res > 1.5) < 0.6)
  expect_true(mean(res > 1.5) > 0.4)
  expect_true(mse(res, dat$y) < 0.1)
})

test_that("Should work without normalization", {
  dat <- rbind(data.frame(y=2, x1=rnorm(500, 1), x2=rnorm(500, 1)),
               data.frame(y=1, x1=rnorm(500,-1), x2=rnorm(500,-1)))
  mod <- LinearizedSVRTrain(X=as.matrix(dat[-1]), Y=dat$y, nump=6, scale=FALSE)
  res <- predict(mod, newdata=as.matrix(dat[-1]))
  expect_true(mean(res > 1.5) < 0.6)
  expect_true(mean(res > 1.5) > 0.4)
  expect_true(mse(res, dat$y) < 0.15)  # Result just isn't as good
})


## An example that goes very badly with clusterY=FALSE, but often
## (not always) does well with clusterY=TRUE:
set.seed(123)
test_that("Clustering with target variable works", {
  dat2 <- rbind(data.frame(y=0, x1=rnorm(500,0), x2=rnorm(500,3)),
                data.frame(y=0, x1=rnorm(500,0), x2=rnorm(500,-3)))
  dat2$y <- 1+(dat2$x1 > 0)
  mod2 <- LinearizedSVRTrain(X=as.matrix(dat2[-1]), Y=dat2$y, nump=2, clusterY=TRUE)
  res2 <- predict(mod2, newdata=as.matrix(dat2[-1]))

  expect_true(mean(res2 > 1.5) < 0.6)
  expect_true(mean(res2 > 1.5) > 0.4)
  expect_true(mse(res2>1.5, dat2$y>1.5) < 0.1)
})
