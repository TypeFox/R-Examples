###################################################
### Chap06Start
###################################################
library(mistat)


###################################################
### SampleFunction
###################################################
set.seed(123)

sample(x=100, size=10)


###################################################
### SampleFromVector
###################################################
X <- 1:100

XSample <- sample(X, size=20, replace=TRUE)


###################################################
### BootMean
###################################################
library(boot)

set.seed(123)

B <- boot(data=X, 
          statistic=function(x, i, n){
            mean(x[i[1:n]])
          },
          R=1000, n=20)

head(B$t, 3)

table(cut(B$t, 12))


###################################################
### FigHistogram1000SampleMeans
###################################################
hist(B$t, main="", xlab="")


###################################################
### RSWORMeanConf
###################################################
X <- 1:100

set.seed(123)

XSmp <- replicate(1000, sample(X, 
                               size=30, 
                               replace=FALSE))

Confint <- function(x, p, n=length(x), N){
  p <- if(p >= 0.5)
    1-((1-p)/2)
  else 
    1-(p/2)
  m <- mean(x)
  z <- qnorm(p=p)/sqrt(n)*(1-((n-1)/(N-1)))^(1/2)
  s <- sd(x)
  res <- m - z*s
  res <- c(res, m + z*s)
  names(res) <- c("lower", "upper")
  return(res)
}

XSmpCnf <- t(apply(XSmp, MARGIN=2, 
                   FUN=Confint, 
                   p=0.95, 
                   N=100))

head(XSmpCnf, 3)

sum(apply(XSmpCnf, MARGIN=1, 
          FUN=function(x, m){
            x[1]< m && x[2] > m
          }, 
          m =50.5))/nrow(XSmpCnf)


###################################################
### FigPredictionSampleMean
###################################################
data(PRED)

set.seed(123)

YPred <- boot(data=PRED$x, 
              statistic=function(x,i){
                mean(x[i[1:100]])* 0.05
              }, 
              R=1000)$t

hist(YPred, main="", 
     xlab="", 
     xlim=c(7, 8))


###################################################
### FigPredictionRatio
###################################################
set.seed(123)

YRatioPred <- boot(data=PRED$x, 
                   statistic=function(x,i){
                     mean(x[i[1:100]])*7.495/148.58
                   }, 
                   R=1000)$t

hist(YRatioPred, main="", 
     xlab="", 
     xlim=c(7, 8))


###################################################
### Chap06End
###################################################
rm(PRED, XSmp, XSmpCnf, YPred, YRatioPred, 
   B, X, XSample, Confint)
detach(package:boot)
detach(package:mistat)
