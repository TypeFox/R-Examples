rliu <-
function(formula,r,R,delt,d,data=NULL,na.action,...)
{
d<-as.matrix(d)
d1<-d[1L]
rliues<-function(formula,r,R,delt,d1,data=NULL,na.action,...)
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
del<-as.matrix(delt)
bb<-xin%*%t(x)%*%y
rb<-bb+solve(s)%*%t(RR)%*%solve(RR%*%solve(s)%*%t(RR))%*%(r-RR%*%bb)
I<-diag(NCOL(x))
fd<-solve(s+I)%*%(s+d1*I)
brs<-fd%*%rb
colnames(brs) <- c("Estimate")
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
dbd<-ev*fd%*%(solve(s)-solve(s)%*%t(RR)%*%solve(RR%*%solve(s)%*%t(RR))%*%RR%*%solve(s))%*%t(fd)
bibet<-(d1-1)*fd%*%solve(s+d1*I)%*%bb+fd%*%t(RR)%*%solve(RR%*%solve(s)%*%t(RR))%*%del
bibets<-bibet%*%t(bibet)
mse<-dbd+bibets
mse1<-sum(diag(mse))
Standard_error<-sqrt(diag(abs(dbd)))
rdel<-matrix(del,nrow(RR))
lenr<-length(RR)
dlpt<-diag(RR%*%xin%*%t(RR))
if (lenr==ncol(RR)) ilpt<-sqrt(solve(abs(dlpt))) else ilpt<-sqrt(solve(diag(abs(dlpt))))
upt<-RR%*%brs
tb<-t(upt)
t_statistic<-((tb-t(rdel))%*%ilpt)/sqrt(ev)
tst<-t(2L*pt(-abs(t_statistic),df<-(NROW(x)-NCOL(x))))
pvalue<-c(tst,rep(NA,(NCOL(x)-NROW(RR))))
mse1<-sum(diag(mse))
mse1<-round(mse1, digits <- 4L)
names(mse1)<-c("MSE")
t_statistic<-c(t_statistic,rep(NA,(NCOL(x)-NROW(RR))))
ans1<-cbind(brs,Standard_error,t_statistic,pvalue)
ans<-round(ans1, digits <- 4L)
anw<-list("*****Restricted Liu Estimator*****"=ans,"*****Mean square error value*****"=mse1)
return(anw)
}
npt<-rliues(formula,r,R,delt,d1,data,na.action)
plotrliu<-function(formula,r,R,delt,d,data=NULL,na.action,...)
{
i<-0
arr<-0
for (i in 1:NROW(d))
{
if (d[i]<0L) d[i]<-0L else d[i]<-d[i]
if (d[i]>1L) d[i]<-1L else d[i]<-d[i]
rlim<-function(formula,r,R,delt,d,data,na.action,...)
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
if (is.matrix(R)) RR<-R else R<-RR
del<-as.matrix(delt)
bb<-xin%*%t(x)%*%y
I<-diag(NCOL(x))
fd<-solve(s+I)%*%(s+d*I)
rb<-bb+solve(s)%*%t(RR)%*%solve(RR%*%solve(s)%*%t(RR))%*%(r-RR%*%bb)
brs<-fd%*%rb
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
dbd<-ev*fd%*%(solve(s)-solve(s)%*%t(RR)%*%solve(RR%*%solve(s)%*%t(RR))%*%RR%*%solve(s))%*%t(fd)
bibet<-(d-1)*fd%*%solve(s+d*I)%*%bb+fd%*%t(RR)%*%solve(RR%*%solve(s)%*%t(RR))%*%del
bibets<-bibet%*%t(bibet)
mse1<-dbd+bibets
mse<-sum(diag(mse1))
return(mse)
}
arr[i]<-rlim(formula,r,R,delt,d[i],data,na.action)
}
MSE<-arr
parameter<-d
pvl<-cbind(parameter,MSE)
colnames(pvl)<-c("Parameter","MSE")
sval<-pvl[order(MSE),]
return(sval)
}
prliu<-plotrliu(formula,r,R,delt,d,data,na.action)
if (nrow(d)>1L) val<-prliu else val<-npt
val
}
