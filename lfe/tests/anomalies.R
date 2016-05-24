library(lfe)
set.seed(42)
options(digits=3, warn=1,lfe.threads=1)
x <- rnorm(1000)
x1 <- rnorm(1000)
x2 <- rnorm(1000) 
z1 <- rnorm(1000)
z2 <- rnorm(1000)
Q <- rnorm(1000) 
f1 <- sample(8,length(x), repl=T)
f2 <- factor(sample(8,length(x), repl=T))
clu <- factor(sample(30,length(x), replace=T))
cluerr <- rnorm(nlevels(clu))[clu]
err <- abs(x)*rnorm(length(x)) + cluerr
y <- x +rnorm(nlevels(clu),sd=0.3)[clu] +  log(f1) + err
f1 <- factor(f1)
dat <- data.frame(y, x, f1=f1, f2, cluster=clu)

# deprecated stuff
summary(felm(y ~x + f1+f2, dat, clustervar='clu'))
summary(felm(y ~ x + G(f1)+G(f2), dat, clustervar=clu))
summary(felm(y ~ x | f1+f2, dat, clustervar=clu))

#anomalies. No variables, etc.
summary(est <- felm(y ~ 1 | f1|0|cluster, dat))
fevcov(est)
summary(felm(y ~ 0 | f1|0|cluster, dat))
summary(felm(y ~ 0 | 0|0|cluster, dat))

summary(felm(y ~ x + x2|f1+factor(f2),dat))
summary(felm(y ~ x + x2+f1|factor(f2),dat))
summary(felm(y ~ x + x2+f1+factor(f2),dat))
summary(lm(y ~ x + x2 + f1 + factor(f2),dat))
summary(felm(y ~ x + x2 + f1 |0|0|0|factor(f2)))
summary(felm(y ~ x + f1+factor(f2) |0|0|0|x2))

# IV
est <- felm(y ~ x | 0 | (x1 | x2 ~ z1 + z2))
for(lh in est$stage1$lhs) print(summary(est$stage1, lhs=lh))
summary(est)
condfstat(est,type=NULL)
summary(est2 <- felm(y ~1 | 0 | (x1 | x2 ~ z1 + z2) | 0 | x))
condfstat(est2, type=NULL)

est0 <- felm( y ~ 1|0|(Q~z1))
condfstat(est0)

# autoload plm:
if(require('plm', quietly=TRUE)) {
  data("EmplUK", package = "plm")
  Em <- pdata.frame(EmplUK)
  detach('package:plm', unload=TRUE)
  print(felm(emp~output+capital + lag(wage,1)|firm, data=Em))
}
