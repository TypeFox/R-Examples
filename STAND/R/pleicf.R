pleicf <-
function(dd,nondet=TRUE,mine=1e-06,maxc=10000,eps=1e-14){
if(nondet){
left<- ifelse(dd[,2]==1, dd[,1],0 )
right<- dd[,1]
}
else { left<-dd[,1] ; right<-dd[,2] }
n<- dim(dd)[1]

L<-ifelse(left==right,left-eps/2,left+eps)
R<-ifelse(left < right,right-eps,right)
n<-length(R)
 o<- icfit(L,R,minerror=mine,maxcount=maxc )
 # ot<<- o
 a<- o$time ; surv<-o$surv ; ple<- 1- surv 
 p<- o$p; p<- p[1:length(p)-1] 
 out<-  cbind(a,ple,surv,p)
 out<- out[ out[,4]!=0,]
 dimnames(out)<-list(NULL,c("a","ple","surv","prob") )
 out<- data.frame(out,n)
out
}

