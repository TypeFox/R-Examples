aul <-
function(formula,d,data=NULL,na.action,...)
{     
d<-as.matrix(d)
d1<-d[1L]
aules<-function(formula,d1,data=NULL,na.action)
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
bb<-xin%*%t(x)%*%y
I<-diag(NCOL(x))
fd<-solve(s+I)%*%(s+d1*I)
bal<-(I-(1-d1)^2*solve(s+I)%*%solve(s+I))%*%bb
colnames(bal) <- c("Estimate")
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
dbd<-ev*(I+(1-d1)*solve(s+I))%*%fd%*%xin%*%fd%*%(I+(1-d1)*solve(s+I))
Standard_error<-sqrt(diag(abs(dbd)))
dbt<-t(bal)
dbd<-ev*(I+(1-d1)*solve(s+I))%*%fd%*%xin%*%fd%*%(I+(1-d1)*solve(s+I))
sdbd_inv<-(sqrt(diag(abs(dbd))))^-1
sdbd_inv_mat<-diag(sdbd_inv)
if (NCOL(dbt)==1L) tbd<-dbt*sdbd_inv else tbd<-dbt%*%sdbd_inv_mat
hggh<-t(tbd)
bibet<--(1-d1)^2*solve(s+I)%*%solve(s+I)%*%bb
bibets<-bibet%*%t(bibet)
mse<-dbd+bibets
mse1<-sum(diag(mse))
mse1<-round(mse1, digits = 4L)
names(mse1)<-c("MSE")
tst<-t(2L*pt(-abs(tbd),df<-(NROW(x)-NCOL(x))))
        colnames(tst) <- c("p_value")
colnames(hggh) <- c("t_statistic")
ans1<-cbind(bal,Standard_error,hggh,tst)
ans<-round(ans1, digits = 4L)
anw<-list("*****Almost Unbiased Liu Estimator*****"=ans,"*****Mean square error value*****"=mse1)
return(anw)
}
npt<-aules(formula,d1,data,na.action)
plotaul<-function(formula,d,data=NULL,na.action,...)
{
i<-0
arr<-0
for (i in 1:NROW(d))
{
if (d[i]<0L) d[i]<-0L else d[i]<-d[i]
if (d[i]>1L) d[i]<-1L else d[i]<-d[i]
aulm<-function(formula,d,data,na.action,...)
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
bb<-xin%*%t(x)%*%y
I<-diag(NCOL(x))
fd<-solve(s+I)%*%(s+d*I)
bal<-(I-(1-d)^2*solve(s+I)%*%solve(s+I))%*%bb
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
dbd<-ev*(I+(1-d)*solve(s+I))%*%fd%*%xin%*%fd%*%(I+(1-d)*solve(s+I))
bibet<--(1-d)^2*solve(s+I)%*%solve(s+I)%*%bb
bibets<-bibet%*%t(bibet)
mse<-dbd+bibets
mse1<-sum(diag(mse))
return(mse1)
}
arr[i]<-aulm(formula,d[i],data,na.action)
}
MSE<-arr
parameter<-d
pvl<-cbind(parameter,MSE)
colnames(pvl)<-c("Parameter","MSE")
sval<-pvl[order(MSE),]
return(sval)
}
paul<-plotaul(formula,d,data,na.action)
if (nrow(d)>1L) val<-paul else val<-npt
val
}
