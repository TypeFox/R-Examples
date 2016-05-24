rid <-
function(formula,k,data=NULL,na.action,...)
{
k<-as.matrix(k)
k1<-k[1L]
rides<-function(formula,k1,data=NULL,na.action)
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
I<-diag(NCOL(x))
bb<-solve(s)%*%t(x)%*%y
bk<-solve(s+k1*I)%*%t(x)%*%y
colnames(bk) <- c("Estimate")
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
dbd<-ev*solve(s+k1*I)%*%s%*%solve(s+k1*I)
Standard_error<-sqrt(diag(abs(dbd)))
dbt<-t(bk)
dbd<-ev*solve(s+k1*I)%*%s%*%solve(s+k1*I)
sdbd_inv<-(sqrt(diag(abs(dbd))))^-1
sdbd_inv_mat<-diag(sdbd_inv)
if (NCOL(dbt)==1L) tbd<-dbt*sdbd_inv else tbd<-dbt%*%sdbd_inv_mat
hggh<-t(tbd)
bibet<--k1*solve(s+k1*I)%*%bb
bibets<-bibet%*%t(bibet)
mse<-dbd+bibets
mse1<-sum(diag(mse))
mse1<-round(mse1, digits <- 4L)
names(mse1)<-c("MSE")
tst<-t(2L*pt(-abs(tbd),df<-(NROW(x)-NCOL(x))))
        colnames(tst) <- c("p_value")
colnames(hggh) <- c("t_statistic")
ans1<-cbind(bk,Standard_error,hggh,tst)
ans<-round(ans1, digits <- 4L)
anw<-list("*****Ordinary Ridge Regression Estimators*****"=ans,"*****Mean square error value*****"=mse1)
return(anw)
}
npt<-rides(formula,k1,data,na.action)
plotrid<-function(formula,k,data=NULL,na.action,...)
{
j<-0
arr<-0
for (j in 1:nrow(k))
{
ridm<-function(formula,k,data,na.action,...)
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
I<-diag(NCOL(x))
bb<-solve(s)%*%t(x)%*%y
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
dbd<-ev*solve(s+k*I)%*%s%*%solve(s+k*I)
bibet<-k*solve(s+k*I)%*%bb
bibets<-bibet%*%t(bibet)
mse<-dbd+bibets
mses<-sum(diag(mse))
return(mses)
}
arr[j]<-ridm(formula,k[j],data,na.action)
}
MSE<-arr
parameter<-k
pvl<-cbind(parameter,MSE)
colnames(pvl)<-c("Parameter","MSE")
sval<-pvl[order(MSE),]
return(sval)
}
prdes<-plotrid(formula,k,data,na.action)
if (nrow(k)>1L) val<-prdes else val<-npt
val
}
