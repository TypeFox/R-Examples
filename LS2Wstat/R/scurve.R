scurve <-function(x,a=1,start=1,end=2){

x<-x*20 - 10

b<-end-start

y<-1/(1+exp(-a*x))

y<-y*b+start

return(y)
}

