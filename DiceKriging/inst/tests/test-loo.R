
x <- c(0, 0.4, 0.6, 0.8, 1)
y <- c(-0.3, 0, -0.8, 0.5, 0.9)
n <- length(x)

theta <- 0.01; sigma <- 3;

list.type <- list("SK", "SK", "UK", "UK")
list.trend.reestim <- list(FALSE, TRUE, FALSE, TRUE)

list.case <- list("", "with known nugget") #, "with noisy observations")
list.nugget <- list(NULL, sigma/10) #, NULL)
list.noise.var <- list(NULL, NULL) #, 1:n)

precision <- 1e-6


for (b in 1:length(list.type)) {
  if (list.trend.reestim[[b]]) {
    trend <- NULL
  } else {
    trend <- c(-1,2)
  }
    
for (a in 1:length(list.case)) {
  
  case <- paste(list.type[[b]], list.case[[a]], "; trend.reestim=", list.trend.reestim[[b]])
  print(case)
  
  m <- km(~x, design=data.frame(x=x), response=data.frame(y=y), 
          covtype="matern5_2", coef.trend=trend, coef.cov=theta, coef.var=sigma^2,
          nugget=list.nugget[[a]], noise.var=list.noise.var[[a]])

  loo.mean <- loo.sd <- rep(NA, n)
  for (i in 1:n){
    mloo <- km(~x, design=data.frame(x=x[-i]), response=data.frame(y=y[-i]), 
               covtype="matern5_2", coef.trend=trend, coef.cov=theta, coef.var=sigma^2,
               nugget=list.nugget[[a]], noise.var=list.noise.var[[a]][-i])
    p <- predict(mloo, newdata=x[i], type=list.type[[b]], checkNames=FALSE)
    loo.mean[i] <- p$mean
    loo.sd[i] <- p$sd
  }

  loo <- leaveOneOut.km(m, type=list.type[[b]], trend.reestim=list.trend.reestim[[b]])
  
  test_that(desc=paste("Check LOO mean for", list.type[[b]], list.case[[a]]), expect_true(max(abs(loo$mean-loo.mean)/loo.mean) < precision))
  test_that(desc=paste("Check LOO sd for", list.type[[b]], list.case[[a]]), expect_true(max(abs(loo$sd-loo.sd)/loo.sd) < precision))
  sum((loo$mean - loo.mean)^2)
  sum((loo$sd - loo.sd)^2)
}
  
}
