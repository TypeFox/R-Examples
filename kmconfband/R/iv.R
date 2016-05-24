iv<-function(sobj){

k<-sum(sobj$n.event>0)

if (k <= 100) q95<-(3.6792+0.5720*log(k)-0.0567*(log(k))^2+0.0027*(log(k))^3)/k
   else       q95<-(3.7752+0.5062*log(k)-0.0417*(log(k))^2+0.0016*(log(k))^3)/k

q95}
