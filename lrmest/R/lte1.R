lte1 <-
function(formula,k,d,press=FALSE,data=NULL,na.action,...)
{
k<-as.matrix(k)
d<-as.matrix(d)
k1<-k[1L]
d1<-d[1L]
ltes1<-function(formula,k1,d1,data=NULL,na.action,...)
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
bb<-solve(t(x)%*%x)%*%t(x)%*%y
bk<-(solve(s+k1*I)-d1*solve(s+k1*I)%*%solve(s))%*%t(x)%*%y
colnames(bk) <- c("Estimate")
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
dbd<-ev*((solve(s+k1*I)-d1*solve(s+k1*I)%*%solve(s))%*%s%*%(solve(s+k1*I)-d1*solve(s)%*%solve(s+k1*I)))
Standard_error<-sqrt(diag(abs(dbd)))
dbt<-t(bk)
dbd<-ev*((solve(s+k1*I)-d1*solve(s+k1*I)%*%solve(s))%*%s%*%(solve(s+k1*I)-d1*solve(s)%*%solve(s+k1*I)))
sdbd_inv<-(sqrt(diag(abs(dbd))))^-1
sdbd_inv_mat<-diag(sdbd_inv)
if (NCOL(dbt)==1L) tbd<-dbt*sdbd_inv else tbd<-dbt%*%sdbd_inv_mat
hggh<-t(tbd)
tst<-t(2L*pt(-abs(tbd),df<-(NROW(x)-NCOL(x))))
        colnames(tst) <- c("p_value")
colnames(hggh) <- c("t_statistic")
bibet<-(solve(s+k1*I)%*%s-d1*solve(s+k1*I)-I)%*%bb
bibets<-bibet%*%t(bibet)
mse<-dbd+bibets
mse1<-sum(diag(mse))
mse1<-round(mse1, digits <- 4L)
names(mse1)<-c("MSE")
ans1<-cbind(bk,Standard_error,hggh,tst)
pr<-0
 i<-1
 m<-1
for (i in 1:nrow(x))
{
               subsum<-0
bb<-c(x[i,]%*%t(x[i,]))
 z<-(solve(s+k1*I-bb)-d1*solve(s+k1*I-bb)%*%solve(s-bb))%*%(t(x)%*%y-x[i,]*y[i])
            
for (m in 1:ncol(x))      
subsum<-subsum+(x[i,m]%*%z[m])         
pr<-pr+(y[i]-subsum)^2
}
pre<-pr
colnames(pre)<-c("PRESS")
ans<-round(ans1, digits = 4L)
pre<-round(pre, digits = 4L)
anw<-list("*****Type (1) Liu Estimator*****"=ans,"*****Mean Square Error value*****"=mse1,"*****Prediction Sum of Squares value*****" = pre)
return(anw)
}
npt<-ltes1(formula,k1,d1,data,na.action)
plotlt1<-function(formula,k,d,data=NULL,na.action,...)
{
j<-0
i<-0
arr<-0
for (j in 1:nrow(k))
{
for (i in 1:nrow(d))
{
lt1m<-function(formula,k,d,data,na.action,...)
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
bb<-solve(t(x)%*%x)%*%t(x)%*%y
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
dbd<-ev*(solve(s+k*I)-d*solve(s+k*I)%*%solve(s))%*%s%*%(solve(s+k*I)-d*solve(s)%*%solve(s+k*I))
bibet<-(solve(s+k*I)%*%s-d*solve(s+k*I)-I)%*%bb
bibets<-bibet%*%t(bibet)
mse<-dbd+bibets
msem<-sum(diag(mse))
return(msem)
}
arr[j*i]<-lt1m(formula,k[j],d[i],data,na.action)
flte1<-file("lte1v.data","a+")
cat(k[j],d[i],arr[j*i],"\n",file=flte1,append=TRUE)
close(flte1)
}
}
mat<-read.table("lte1v.data")
unlink("lte1v.data")
rmat<-matrix(mat[,3L],nrow=NROW(d),dimnames=list(c(paste0("d=",d)),c(paste0("k=",k))))
return(rmat)
}
plte1<-plotlt1(formula,k,d,data,na.action)
plotprlt1<-function(formula,k,d,data=NULL,na.action,...)
{
j<-0
i<-0
arr<-0
for (j in 1:nrow(k))
{
for (i in 1:nrow(d))
{
pre<-function(formula,k,d,data,na.action,...)
{
cal<-match.call(expand.dots=FALSE)
mat<-match(c("formula","data","na.action"), names(cal))
cal<-cal[c(1,mat)]
cal[[1]]<-as.name("model.frame") 
cal<-eval(cal)
y<-model.response(cal)         
md<- attr(cal, "terms")
x<-model.matrix(md,cal,contrasts)
s<-t(x)%*%x
I<-diag(NCOL(x))
pr<-0
 i<-1
 m<-1
for (i in 1:nrow(x))
{
               subsum<-0
bb<-c(x[i,]%*%t(x[i,]))
 z<-(solve(s+k*I-bb)-d*solve(s+k*I-bb)%*%solve(s-bb))%*%(t(x)%*%y-x[i,]*y[i])
            
for (m in 1:ncol(x))      
subsum<-subsum+(x[i,m]%*%z[m])         
pr<-pr+(y[i]-subsum)^2
}
pre<-t(pr)
return(pre)
}
arr[j*i]<-pre(formula,k[j],d[i],data,na.action)
fprlt1<-file("prltv1.data","a+")
cat(k[j],d[i],arr[j*i],"\n",file=fprlt1,append=TRUE)
close(fprlt1)
}
}
pmat<-read.table("prltv1.data")
unlink("prltv1.data")
rprmat<-matrix(pmat[,3L],nrow=NROW(d),dimnames=list(c(paste0("d=",d)),c(paste0("k=",k))))
return(rprmat)
}
prlt1<-plotprlt1(formula,k,d,data,na.action)
if (!press) prmse<-plte1 else prmse<-prlt1
if (nrow(k)>1L | nrow(k)>1L) val<-prmse else val<-npt
val
}
