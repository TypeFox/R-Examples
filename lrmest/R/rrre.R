rrre <-
function(formula,r,R,dpn,delt,k,data=NULL,na.action,...)
{
k<-as.matrix(k)
k1<-k[1L]
rrres<-function(formula,r,R,dpn,delt,k1,data=NULL,na.action,...)
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
del<-delt
de1<-as.matrix(delt)
I<-diag(NCOL(x))
wk<-solve(I+k1*s)
bb<-xin%*%t(x)%*%y
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
bk<-wk%*%(bb+solve(s)%*%t(RR)%*%solve(RR%*%solve(s)%*%t(RR))%*%(r-RR%*%bb))
colnames(bk) <- c("Estimate")
dbd<-ev*wk%*%(solve(s)-solve(s)%*%t(RR)%*%(RR%*%solve(s)%*%t(RR))%*%RR%*%solve(s))%*%t(wk)
Standard_error<-sqrt(diag(abs(dbd)))
rdel<-matrix(delt,NROW(RR))
lenr<-length(RR)
dlpt<-diag(RR%*%xin%*%t(RR))
if (lenr==ncol(RR)) ilpt<-sqrt(solve(abs(dlpt))) else ilpt<-sqrt(solve(diag(abs(dlpt))))
upt<-RR%*%bk
tb<-t(upt)
t_statistic<-((tb-t(rdel))%*%ilpt)/sqrt(ev)
tst<-t(2L*pt(-abs(t_statistic),df<-(NROW(x)-NCOL(x))))
pvalue<-c(tst,rep(NA,(NCOL(x)-NROW(RR))))
dbd<-ev*wk%*%(solve(s)-solve(s)%*%t(RR)%*%(RR%*%solve(s)%*%t(RR))%*%RR%*%solve(s))%*%t(wk)
bibet<-wk%*%solve(s)%*%t(RR)%*%solve(RR%*%solve(s)%*%t(RR))%*%del-k1*solve(s+k1*I)%*%bb
bibets<-bibet%*%t(bibet)
mse<-dbd+bibets
mse1<-sum(diag(mse))
mse1<-round(mse1, digits <- 4L)
names(mse1)<-c("MSE")
t_statistic<-c(t_statistic,rep(NA,(NCOL(x)-NROW(RR))))
ans1<-cbind(bk,Standard_error,t_statistic,pvalue)
ans<-round(ans1, digits <- 4L)
anw<-list("*****Restricted Ridge Regression Estimator*****"=ans,"*****Mean square error value*****"=mse1)
return(anw)
}
npt<-rrres(formula,r,R,dpn,delt,k1,data,na.action)
plotrrr<-function(formula,r,R,dpn,delt,k,data=NULL,na.action,...)
{
j<-0
arr<-0
for (j in 1:nrow(k))
{
rrrem<-function(formula,r,R,dpn,delt,k,data,na.action,...)
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
del<-delt
de1<-as.matrix(delt)
I<-diag(NCOL(x))
bb<-xin%*%t(x)%*%y
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
wk<-solve(I+k*s)
dbd<-ev*wk%*%(solve(s)-solve(s)%*%t(RR)%*%(RR%*%solve(s)%*%t(RR))%*%RR%*%solve(s))%*%t(wk)
bibet<-wk%*%solve(s)%*%t(RR)%*%solve(RR%*%solve(s)%*%t(RR))%*%del-k*solve(s+k*I)%*%bb
bibets<-bibet%*%t(bibet)
mse<-dbd+bibets
mse1<-sum(diag(mse))
return(mse1)
}
arr[j]<-rrrem(formula,r,R,dpn,delt,k[j],data,na.action)
}
MSE<-arr
Parameter<-k
pvl<-cbind(Parameter,MSE)
colnames(pvl)<-c("Parameter","MSE")
sval<-pvl[order(MSE),]
return(sval)
}
prrre<-plotrrr(formula,r,R,dpn,delt,k,data,na.action)
if (nrow(k)>1L) val<-prrre else val<-npt
val
}
