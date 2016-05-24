epan<-function(x)
{
d<-length(x)
val<-1
for (i in 1:d){
   val<-val*3*(1-x[i]^2)/2
}
return(val)
}
