`skewness` <-
function(x) {
x <- na.omit(x)
n<-length(x)
s<-sqrt(var(x))
suma<-sum((x-mean(x))^3)/s^3
k <- n*suma/((n-1)*(n-2))
return(k)
}

