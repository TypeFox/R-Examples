SES <-
function(ts,alpha=0.5,s0=NULL){
if(is.null(s0)) s0<-mean(ts)
Q<-rep(0,length(ts))
for(t in 1:length(ts)){
if(t==1)Q[t]<-s0
if(t==2)Q[t]<-alpha*ts[1]+(1-alpha)*s0
if(t==3)Q[t]<-alpha*ts[2]+(1-alpha)*ts[1]+(1-alpha)*s0
if(t>3){
Q[t]=Q[t]+alpha*ts[t-1]+(1-alpha)*ts[t-2]
for(i in 2:(t-2)) Q[t]=Q[t]+alpha*((1-alpha)^i)*ts[t-i-1]
Q[t]=Q[t]+(1-alpha)*s0
}#if t>3
}#for(t in 1:length(ts))
Q
}
