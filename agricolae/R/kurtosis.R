`kurtosis` <-
function(x) {
x<-na.omit(x)
n<-length(x)
suma<-sum((x-mean(x))^4)/(var(x))^2
k <- n*(n+1)*suma/((n-1)*(n-2)*(n-3)) - 3*(n-1)^2/((n-2)*(n-3))
return(k)
}

