data <- c(18.0,6.3,7.5,8.1,3.1,0.8,2.4,3.5,9.5,39.7,
          3.4,14.6,5.1,6.8,2.6,8.0,8.5,3.7,21.2,3.1,
          10.2, 8.3,6.4,3.0,5.7,5.6,7.4,3.9,9.1,4.0)
n <- length(data)
theta.hat <- 1/mean(data); theta.hat
cutpts <- c(0,2,4,7,12,Inf)
bin.data <- cut(data, cutpts)
p <- diff(pexp(cutpts,theta.hat))
e <- n * p
o <- table(bin.data)
print(cbind(o,e))
lrt  <- 2 * sum(o * log(o/e)); lrt
pearson <- sum( (o-e)^2/e ); pearson
1-pchisq(lrt, df=2)               # df = (4 - 1) - 1 [anti-conservative]
1-pchisq(pearson,df=2)
1-pchisq(lrt, df=3)               # df = 4 - 1       [conservative]
1-pchisq(pearson,df=3)
