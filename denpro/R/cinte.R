cinte<-function(values,volumes,parents){
#Calculates the integral of a piecewise continuous function.
#
len<-length(values)
int<-0
for (i in len:1){
  par<-parents[i]
  if (par==0) valpar<-0 
  else valpar<-values[par]
  int<-int+volumes[i]*(values[i]-valpar)
}
return(int)
}

