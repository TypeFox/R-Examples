rgbinom <-
function(N,size,prob){
z<-runif(N,0,1)
x<-qgbinom(z,size,prob)
return(x)
}
