noe.compute.pv<-function(tn,tc){
# Compute the probability vector assuming a Uniform[0,1] cdf
tp<-rep(0,2*tn+1)
tp[1]<-punif(tc[2])
m<-2
while (m<=(2*tn)){
       tp[m]<-punif(tc[m+1])-punif(tc[m])
       m<-m+1}
tp[2*tn+1]<-1-punif(tc[2*tn+1])
tp}
