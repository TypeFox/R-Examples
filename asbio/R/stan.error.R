stan.error<-function(x){
n<-nrow(as.matrix(x))
sqrt(((n*(n-1))^-1)*sum((x-mean(x))^2))}

stan.error.sq<-function(x){
n<-nrow(as.matrix(x))
((n*(n-1))^-1)*sum((x-mean(x))^2)}