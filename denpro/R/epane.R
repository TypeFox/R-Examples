epane<-function(x,h){
#
d<-length(x)
val<-1
for (i in 1:d){
   val<-val*3*(1-(x[i]/h)^2)/2
}
return(val)
}
