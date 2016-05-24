######################################################
# Section 4.2 Normal Data with Both Parameters Unknown
######################################################

library(LearnBayes)

 data(marathontimes)
 attach(marathontimes)
 d = mycontour(normchi2post, c(220, 330, 500, 9000), time,
    xlab="mean",ylab="variance")

 S = sum((time - mean(time))^2) 
 n = length(time)
 sigma2 = S/rchisq(1000, n - 1)
 mu = rnorm(1000, mean = mean(time), sd = sqrt(sigma2)/sqrt(n))

 points(mu, sigma2)

 quantile(mu, c(0.025, 0.975))

 quantile(sqrt(sigma2), c(0.025, 0.975))
