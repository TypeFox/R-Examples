etais<-function(x,y){
#laskee euklid etais nelion vektorien x ja y valilla
#
pit<-length(x)
vast<-0
i<-1
while (i<=pit){
  vast<-vast+(x[i]-y[i])^2
  i<-i+1
}
return(vast)
}
