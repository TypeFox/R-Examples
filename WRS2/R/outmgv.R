outmgv<-function(x,y=NULL,plotit=TRUE,outfun=outbox,se=TRUE,op=1,
cov.fun=rmba,xlab="X",ylab="Y",SEED=TRUE,STAND=FALSE,...){
#
# Check for outliers using mgv method
#
# NOTE: if columns of the input matrix are reordered, this can
# have an effect on the results due to rounding error when calling
# the R function eigen.
#
#  (Argument STAND is included simply to avoid programming issues when outmgv is called by other functions.)
#
if(is.null(y[1]))m<-x
if(!is.null(y[1]))m<-cbind(x,y)
m=elimna(m)
m=as.matrix(m)
nv=nrow(m)
temp<-mgvar(m,se=se,op=op,cov.fun=cov.fun,SEED=SEED)
#if(fast)temp<-mgvdep.for(m,se=se)$distance
temp[is.na(temp)]<-0
if(ncol(m)==1){
temp2=outpro(m)
nout=temp2$n.out
keep=temp2$keep
temp2=temp2$out.id
}
if(ncol(m)>1){
if(ncol(m)==2)temp2<-outfun(temp,...)$out.id
if(ncol(m)>2)temp2<-outbox(temp,mbox=TRUE,gval=sqrt(qchisq(.975,ncol(m))))$out.id
vec<-c(1:nrow(m))
flag<-rep(T,nrow(m))
flag[temp2]<-F
vec<-vec[flag]
vals<-c(1:nrow(m))
keep<-vals[flag]
if(plotit && ncol(m)==2){
x<-m[,1]
y<-m[,2]
plot(x,y,type="n",xlab=xlab,ylab=ylab)
flag<-rep(T,length(y))
flag[temp2]<-F
points(x[flag],y[flag],pch="*")
points(x[temp2],y[temp2],pch="o")
}
nout=0
if(!is.na(temp2[1]))nout=length(temp2)
}
list(n=nv,n.out=nout,out.id=temp2,keep=keep)
}
