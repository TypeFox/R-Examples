cumul.sum <-
function(x){
y<-rep(0,length(x))
#cumulatively sums a vector x and returns the result
for(i in 1:length(x)){
y[i]<-sum(x[1:i])
}
return(y)
}

