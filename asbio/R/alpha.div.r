alpha.div<-function(x,index){
if(index=="simp")alpha<-Simp.index(x)
if(index=="inv.simp")alpha<-Simp.index(x,inv=TRUE)
if(index=="shan")alpha<-SW.index(x)
alpha
}

Simp.index<-function(x,inv=FALSE){
if(ncol(as.matrix(x))==1){
p.i<-x/sum(x)
D<-1-sum(p.i^2)}
if(ncol(as.matrix(x))>1){
p.i<-apply(x,1,function(x){x/sum(x)})
if(inv==FALSE)D<-1-apply(p.i^2,2,sum)
if(inv==TRUE)D<-1/apply(p.i^2,2,sum)
}
D
}

SW.index<-function(x){
if(ncol(as.matrix(x))==1){
p.i<-x/sum(x)
p.i.new<-p.i[p.i!=0]
h.prime<--1*sum(log(p.i.new)*p.i.new)}
if(ncol(as.matrix(x))>1){
p.i<-apply(x,1,function(x){x/sum(x)})
h<-apply(p.i,1,function(x){log(x)*x})
h.prime<- -1*apply(h,1,function(x){sum(x[!is.na(x)])})}
h.prime
}



