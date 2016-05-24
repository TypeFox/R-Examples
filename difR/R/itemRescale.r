# Rescaling of the item parameters
# mR: reference group, mF: focal group
# Item parameters rescaled to reference group metric

itemRescale<-function(mR,mF,items=1:nrow(mR)){
if (ncol(mR)==2) {
res<-cbind(mF[,1]+mean(mR[items,1])-mean(mF[items,1]),mF[,2])
colnames(res)<-c("new.b","se(b)")
}
else{
my<-mean(mR[items,2])
sy<-sd(mR[items,2])
mx<-mean(mF[items,2])
sx<-sd(mF[items,2])
A<-sy/sx
B<-my-mx*sy/sx
new.b<-A*mF[,2]+B
new.a<-mF[,1]/A
if (ncol(mF)==9){
res<-cbind(new.a,new.b,mF[,3],mF[,4]/A,mF[,5]*A,mF[,6:7],mF[,8]/A,mF[,9]*A)
colnames(res)<-c("new.a","new.b","c","new.se(a)","new.se(b)","se(c)","cov(a,b)","new.cov(a,c)","new.cov(b,c)")
}
else{
res<-cbind(new.a,new.b,mF[,3]/A,mF[,4]*A,mF[,5])
if (ncol(mF)==5) colnames(res)<-c("new.a","new.b","new.se(a)","new.se(b)","cov(a,b)")
else{
res<-cbind(res,mF[,6])
colnames(res)<-c("new.a","new.b","new.se(a)","new.se(b)","cov(a,b)","c")
}
}
}
return(res)}