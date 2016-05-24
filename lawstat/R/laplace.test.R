`laplace.test` <-
function(y)
{
### load the VGAM package ###

#library(VGAM)

### omit the NAs ###

 y<-y[!is.na(y)]

### sort NAs ###

 y<-sort(y)
 n<-length(y)
 a<-median(y)
 b<-mean(abs(y-a))
 z <- plaplace((y - a)/b)

### Anderson-Darling statistic ###

 A2 <- -mean((2 * seq(1:n) - 1) * (log(z) + log(1 - rev(z))))-n

### Cramer-von Mises statistic ###

 W2<-sum((z-(2*seq(1:n)-1)/(2*n))^2)+1/(12*n)

### Watson statistic ###

 U2<-W2-n*(mean(z)-0.5)^2

### Kolmogorov statistics (D and V) ###

 D<-sqrt(n)*max(max(seq(1:n)/n-z), max(z-(seq(1:n)-1)/n))
 V<-sqrt(n)*(max(seq(1:n)/n-z)+max(z-(seq(1:n)-1)/n))

### display output ###

 list(A2=A2, W2=W2, U2=U2, D=D, V=V)
}

