rcn<-function(n,eps,sigmac) rnorm(n,0,sample(c(1,sigmac),n,replace=TRUE,prob=c(1-eps,eps)))

