H.mean<-function(x){
x<-na.omit(x)
n<-nrow(as.matrix(x))
((1/n)*sum(1/x))^-1}