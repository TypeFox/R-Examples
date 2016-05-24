srliu <-
function(formula,r,R,dpn,delt,d,data=NULL,na.action,...)
{
d<-as.matrix(d)
d1<-d[1L]
srliues<-function(formula,r,R,dpn,delt,d1,data=NULL,na.action,...)
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
de1t<-as.matrix(delt)
del<-delt
RC<-matrix(R,NCOL(x))
RR<-t(RC)
I<-diag(NCOL(x))
if (is.matrix(R)) RR<-R else RR<-RR
if (length(dpn)==1L) shi<-dpn
else if (is.matrix(dpn)) shi<-dpn else shi<-diag(dpn)
bb<-xin%*%t(x)%*%y
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
w1<-solve(s/ev+t(RR)%*%solve(shi)%*%RR)
w2<-(t(x)%*%y)/ev+t(RR)%*%solve(shi)%*%r
bm<-w1%*%w2
fd<-solve(s+I)%*%(s+d1*I)
srl<-fd%*%bm
colnames(srl) <- c("Estimate")
g11<-solve(s)%*%t(RR)%*%solve((shi/ev)+RR%*%solve(s)%*%t(RR))%*%RR%*%solve(s)
dbd<-ev*fd%*%(solve(s)-g11)%*%t(fd)
Standard_error<-sqrt(diag(abs(dbd)))
dbd<-ev*fd%*%(solve(s)-g11)%*%t(fd)
ta<-(d1-1)*solve(s+d1*I)%*%bb+g11%*%s%*%ginv(RR)%*%del
mse1<-ev*fd%*%(solve(s)-g11)%*%t(fd)+fd%*%ta%*%t(ta)%*%t(fd)
mse<-sum(diag(mse1))
rdel<-matrix(delt,NROW(RR))
lenr<-length(RR)
dlpt<-diag(RR%*%xin%*%t(RR))
if (lenr==ncol(RR)) ilpt<-sqrt(solve(abs(dlpt))) else ilpt<-sqrt(solve(diag(abs(dlpt))))
upt<-RR%*%srl
tb<-t(upt)
t_statistic<-((tb-t(rdel))%*%ilpt)/sqrt(ev)
tst<-t(2L*pt(-abs(t_statistic),df<-(NROW(x)-NCOL(x))))
pvalue<-c(tst,rep(NA,(NCOL(x)-NROW(RR))))
mse<-round(mse, digits <- 4L)
names(mse)<-c("MSE")
t_statistic<-c(t_statistic,rep(NA,(NCOL(x)-NROW(RR))))
ans1<-cbind(srl,Standard_error,t_statistic,pvalue)
ans<-round(ans1, digits <- 4L)
anw<-list("*****Stochastic Restricted Liu Estimator*****"=ans,"*****Mean square error value*****"=mse)
return(anw)
}
npt<-srliues(formula,r,R,dpn,delt,d1,data,na.action)
plotsrliu<-function(formula,r,R,dpn,delt,d,data=NULL,na.action,...)
{
i<-0
arr<-0
for (i in 1:NROW(d))
{
if (d[i]<0L) d[i]<-0L else d[i]<-d[i]
if (d[i]>1L) d[i]<-1L else d[i]<-d[i]
srlim<-function(formula,r,R,dpn,del,d,data,na.action,...)
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
de1t<-as.matrix(delt)
del<-delt
RC<-matrix(R,NCOL(x))
RR<-t(RC)
I<-diag(NCOL(x))
if (is.matrix(R)) RR<-R else RR<-RR
if (length(dpn)==1L) shi<-dpn
else if (is.matrix(dpn)) shi<-dpn else shi<-diag(dpn)
bb<-xin%*%t(x)%*%y
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
w1<-solve(s/ev+t(RR)%*%solve(shi)%*%RR)
fd<-solve(s+I)%*%(s+d*I)
g11<-solve(s)%*%t(RR)%*%solve((shi/ev)+RR%*%solve(s)%*%t(RR))%*%RR%*%solve(s)
dbd<-ev*fd%*%(solve(s)-g11)%*%t(fd)
ta<-(d-1)*solve(s+d*I)%*%bb+g11%*%s%*%ginv(RR)%*%del
mse1<-ev*fd%*%(solve(s)-g11)%*%t(fd)+fd%*%ta%*%t(ta)%*%t(fd)
mse<-sum(diag(mse1))
return(mse)
}
arr[i]<-srlim(formula,r,R,dpn,delt,d[i],data,na.action)
}
MSE<-arr
parameter<-d
pvl<-cbind(parameter,MSE)
colnames(pvl)<-c("Parameter","MSE")
sval<-pvl[order(MSE),]
return(sval)
}
psrle<-plotsrliu(formula,r,R,dpn,delt,d,data,na.action)
if (nrow(d)>1L) val<-psrle else val<-npt
val
}
