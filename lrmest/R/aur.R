aur <-
function(formula,k,data=NULL,na.action,...)
{
k<-as.matrix(k)
k1<-k[1L]
aures<-function(formula,k1,data=NULL,na.action,...)
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
tk<-solve(s+k1*I)%*%s
bar<-(I-k1^2*solve(s+k1*I)%*%solve(s+k1*I))%*%bb
colnames(bar)<-c("Estimate")
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
dbb<-ev*(I+k1*solve(s+k1*I))%*%tk%*%xin%*%tk%*%(I+k1*solve(s+k1*I))
bibet<--k1^2*solve(s+k1*I)%*%solve(s+k1*I)%*%bb
bibets<-bibet%*%t(bibet)
mse<-dbb+bibets
mse1<-sum(diag(mse))
mse1<-round(mse1, digits <- 4L)
names(mse1)<-c("MSE")
Standard_error<-sqrt(diag(abs(dbb)))
dbt<-t(bar)
sdbd_inv<-(sqrt(diag(abs(dbb))))^-1
sdbd_inv_mat<-diag(sdbd_inv)
if (NCOL(dbt)==1L) tbb<-dbt*sdbd_inv else tbb<-dbt%*%sdbd_inv_mat
hggh<-t(tbb)
tst<-t(2L*pt(-abs(tbb),df<-(NROW(x)-NCOL(x))))
        colnames(tst) <- c("p_value")
colnames(hggh) <- c("t_statistic")     
ans1<-cbind(bar,Standard_error,hggh,tst)
ans<-round(ans1, digits = 4L)
adw<-list("*****Almost Unbiased Ridge Estimator******"=ans,"*****Mean square error value*****"=mse1)
return(adw)
}
npt<-aures(formula,k1,data,na.action)
plotaur<-function(formula,k,data=NULL,na.action,...)
{
j<-0
arr<-0
for (j in 1:nrow(k))
{
aurm<-function(formula,k,data,na.action,...)
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
tk<-solve(s+k*I)%*%s
bar<-(I-k^2*solve(s+k*I)%*%solve(s+k*I))%*%bb
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
dbb<-ev*(I+k*solve(s+k*I))%*%tk%*%xin%*%tk%*%(I+k*solve(s+k*I))
bibet<--k^2*solve(s+k*I)%*%solve(s+k*I)%*%bb
bibets<-bibet%*%t(bibet)
mse<-dbb+bibets
mse1<-sum(diag(mse))
return(mse1)
}
arr[j]<-aurm(formula,k[j],data,na.action)

}
MSE<-arr
parameter<-k
pvl<-cbind(parameter,MSE)
colnames(pvl)<-c("Parameter","MSE")
sval<-pvl[order(MSE),]
return(sval)
}
paur<-plotaur(formula,k,data,na.action)
if (nrow(k)>1L) val<-paur else val<-npt
val
}
