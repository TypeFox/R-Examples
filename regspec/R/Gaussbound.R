Gaussbound<-function(mode,tau,prob){
q<-1-prob
if(q<1/3){k<-tau*2/3/q^0.5}
if(q>1/3){k<-tau*3^0.5*(1-q)}
mode+c(-1,1)*k
}