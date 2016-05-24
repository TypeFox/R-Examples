"fun.beta1" <-
function(a, b){
index<-as.logical((a>0)*(b>0))

result<-gamma(a[index])*gamma(b[index])/gamma(a[index] + b[index])
result[!index]<-NA
return(result)
}

