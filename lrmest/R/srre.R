srre <-
function(formula,r,R,dpn,delt,k,data=NULL,na.action,...)
{
k<-as.matrix(k)
k1<-k[1L]
srres<-function(formula,r,R,dpn,delt,k1,data=NULL,na.action,...)
{
cal<-match.call(expand.dots=FALSE)
mat<-match(c("formula","data","na.action"), names(cal))
cal<-cal[c(1L,mat)]
cal[[1L]]<-as.name("model.frame") 
cal<-eval(cal)
y<-model.response(cal)         
md<- attr(cal, "terms")
x<-model.matrix(md,cal,contrasts)
s<-t(x)%*%x
xin<-solve(s)
r<-as.matrix(r)
RC<-matrix(R,NCOL(x))
RR<-t(RC)
if (is.matrix(R)) RR<-R else RR<-RR
if (length(dpn)==1L) shi<-dpn
else if (is.matrix(dpn)) shi<-dpn else shi<-diag(dpn)
de1<-as.matrix(delt)
I<-diag(NCOL(x))
bb<-xin%*%t(x)%*%y
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
w1<-solve(s/ev+t(RR)%*%solve(shi)%*%RR)
w2<-(t(x)%*%y)/ev+t(RR)%*%solve(shi)%*%r
bm<-w1%*%w2
tk<-solve(s+k1*I)%*%s
bsrr<-tk%*%bm
colnames(bsrr) <- c("Estimate")
dbd<-ev*tk%*%solve(s+t(RR)%*%solve(shi/ev)%*%RR)%*%tk
Standard_error<-sqrt(diag(abs(dbd)))
rdel<-matrix(delt,NROW(RR))
lenr<-length(RR)
dlpt<-diag(RR%*%xin%*%t(RR))
if (lenr==ncol(RR)) ilpt<-sqrt(solve(abs(dlpt))) else ilpt<-sqrt(solve(diag(abs(dlpt))))
upt<-RR%*%bsrr
tb<-t(upt)
t_statistic<-((tb-t(rdel))%*%ilpt)/sqrt(ev)
tst<-t(2L*pt(-abs(t_statistic),df<-(NROW(x)-NCOL(x))))
pvalue<-c(tst,rep(NA,(NCOL(x)-NROW(RR))))
dbd<-ev*tk%*%solve(s+t(RR)%*%solve(shi/ev)%*%RR)%*%tk
bibet<-k1*solve(s+k1*I)%*%bb
bibets<-bibet%*%t(bibet)
mse<-dbd+bibets
mse1<-sum(diag(mse))
mse1<-round(mse1, digits <- 4L)
names(mse1)<-c("MSE")
t_statistic<-c(t_statistic,rep(NA,(NCOL(x)-NROW(RR))))
ans1<-cbind(bsrr,Standard_error,t_statistic,pvalue)
ans<-round(ans1, digits <- 4L)
anw<-list("*****Stochastic Restricted Ridge Estimator*****"=ans,"*****Mean square error value*****"=mse1)
return(anw)
}
npt<-srres(formula,r,R,dpn,delt,k1,data,na.action)
plotsrr<-function(formula,r,R,dpn,delt,k,data=NULL,na.action,...)
{
j<-0
arr<-0
for (j in 1:nrow(k))
{
srrem<-function(formula,r,R,dpn,delt,k,data,na.action,...)
{
cal<-match.call(expand.dots=FALSE)
mat<-match(c("formula","data","na.action"), names(cal))
cal<-cal[c(1L,mat)]
cal[[1L]]<-as.name("model.frame") 
cal<-eval(cal)
y<-model.response(cal)         
md<- attr(cal, "terms")
x<-model.matrix(md,cal,contrasts)
s<-t(x)%*%x
xin<-solve(s)
r<-as.matrix(r)
RC<-matrix(R,NCOL(x))
RR<-t(RC)
if (is.matrix(R)) RR<-R else RR<-RR
if (length(dpn)==1L) shi<-dpn
else if (is.matrix(dpn)) shi<-dpn else shi<-diag(dpn)
de1<-as.matrix(delt)
I<-diag(NCOL(x))
bb<-xin%*%t(x)%*%y
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
tk<-solve(s+k*I)%*%s
dbd<-ev*tk%*%solve(s+t(RR)%*%solve(shi/ev)%*%RR)%*%tk
bibet<-k*solve(s+k*I)%*%bb
bibets<-bibet%*%t(bibet)
mse<-dbd+bibets
mse1<-sum(diag(mse))
return(mse1)
}
arr[j]<-srrem(formula,r,R,dpn,delt,k[j],data,na.action)
}
MSE<-arr
Parameter<-k
pvl<-cbind(Parameter,MSE)
colnames(pvl)<-c("Parameter","MSE")
sval<-pvl[order(MSE),]
return(sval)
}
psrre<-plotsrr(formula,r,R,dpn,delt,k,data,na.action)
if (nrow(k)>1L) val<-psrre else val<-npt
val
}
