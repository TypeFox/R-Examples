dgbinom <-
function(x,size,prob,log=FALSE){
n<-sum(size)
y<-round(x)
if (any(is.na(x) | (x < 0)) || max(abs(x - y)) > 1e-07) 
        stop("'x' must be nonnegative and integer")
z<-round(n)
 if (length(size) != length(prob) || is.na(n) || (n < 1) || abs(n - 
            z) > 1e-07 || (x > z)) 
            stop("'size' must contain positive integers and their sum must be  >= 'x'")
if (any(is.na(prob)| (prob<0) | (prob>1) ))
stop("'prob' must contain numbers between 0 und 1")
x<-x+1
theta<-c(rep(prob,size))
xi<-c(1-theta[1],theta[1])
if(n>1){
for(i in c(2:n)){
P<-matrix(c(xi,0,0,xi),i+1,2)
xi<-P%*%c(1-theta[i],theta[i])}}
if (log==TRUE){
xi<-log(xi)}
xi[x]
}
