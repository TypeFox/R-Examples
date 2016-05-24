alte1 <-
function(formula,k,d,aa,press=FALSE,data=NULL,na.action,...)
{
k<-as.matrix(k)
d<-as.matrix(d)
k1<-k[1L]
d1<-d[1L]
if (length(aa)==1L) A<-as.matrix(aa) else A<-diag(aa)
altes<-function(formula,k1,d1,aa,data=NULL,na.action,...)
{
cal<-match.call(expand.dots=FALSE)
mat<-match(c("formula","data","na.action"), names(cal))
cal<-cal[c(1L,mat)]
cal[[1L]]<-as.name("model.frame") 
cal<-eval(cal)
y<-model.response(cal,"numeric")
md<- attr(cal, "terms")
x<-model.matrix(md,cal,contrasts)
s<-t(x)%*%x
I<-diag(NCOL(x))
bb1<-solve(s)%*%t(x)%*%y
bk1<-(solve(s+k1*I)-d1*solve(s+k1*I)%*%solve(s))%*%t(x)%*%y
bk<-A%*%bk1
colnames(bk) <- c("Estimate")
ev<-(t(y)%*%y-t(bb1)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
dbd<-ev*A%*%(solve(s+k1*I)-d1*solve(s+k1*I)%*%solve(s))%*%s%*%(solve(s+k1*I)-d1*solve(s)%*%solve(s+k1*I))%*%t(A)
Standard_error<-sqrt(diag(abs(dbd)))
dbt<-t(bk)
dbd<-ev*A%*%(solve(s+k1*I)-d1*solve(s+k1*I)%*%solve(s))%*%s%*%(solve(s+k1*I)-d1*solve(s)%*%solve(s+k1*I))%*%t(A)
sdbd_inv<-(sqrt(diag(abs(dbd))))^-1
sdbd_inv_mat<-diag(sdbd_inv)
if (NCOL(dbt)==1L) tbd<-dbt*sdbd_inv else tbd<-dbt%*%sdbd_inv_mat
hggh<-t(tbd)
tst<-t(2L*pt(-abs(tbd),df<-(NROW(x)-NCOL(x))))
        colnames(tst) <- c("p_value")
colnames(hggh) <- c("t_statistic")
bibet<-A%*%(solve(s+k1*I)%*%s-d1*solve(s+k1*I)-I)%*%bb1
bibets<-bibet%*%t(bibet)
mse<-dbd+bibets
mse1<-sum(diag(mse))
mse1<-round(mse1, digits <- 4L)
names(mse1)<-c("MSE")
ab<-matrix(aa,NCOL(x))
a<-t(ab)
pr<-0
 i<-1
 m<-1
for (i in 1:nrow(x))
{
               subsum<-0
bb<-c(x[i,]%*%t(x[i,]))
 z<-(solve(s+k1*I-bb)-d1*solve(s+k1*I-bb)%*%solve(s-bb))%*%(t(x)%*%y-x[i,]*y[i])
            
for (m in 1:ncol(x))      
subsum<-subsum+(a[m]%*%x[i,m]%*%z[m])         
pr<-pr+(y[i]-subsum)^2
}
pre<-t(pr)
colnames(pre)<-c("PRESS")
ans1<-cbind(bk,Standard_error,hggh,tst)
ans1<-round(ans1, digits = 4L)
pre<-round(pre, digits = 4L)
rownames(ans1)<-rownames(bb1)
anw<-list("*****Type (1) Adjusted Liu Estimator*****"=ans1,"****Mean Square Error value*****"=mse1,"*****Prediction Sum of Squares value*****"=pre)
return(anw)
}
npt<-altes(formula,k1,d1,aa,data,na.action)
plotalt1<-function(formula,k,d,aa,data=NULL,na.action,...)
{
i<-0
j<-0
arr<-0
for (j in 1:nrow(k))
{
for (i in 1:nrow(d))
{
altem1<-function(formula,k,d,aa,data,na.action,...)
{
cal<-match.call(expand.dots=FALSE)
mat<-match(c("formula","data","na.action"), names(cal))
cal<-cal[c(1L,mat)]
cal[[1L]]<-as.name("model.frame") 
cal<-eval(cal)
y<-model.response(cal,"numeric")
md<- attr(cal, "terms")
x<-model.matrix(md,cal,contrasts)
s<-t(x)%*%x
I<-diag(NCOL(x))
bb<-solve(t(x)%*%x)%*%t(x)%*%y
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
dbd<-ev*A%*%(solve(s+k*I)-d*solve(s+k*I)%*%solve(s))%*%s%*%(solve(s+k*I)-d*solve(s)%*%solve(s+k*I))%*%t(A)
bibet<-A%*%(solve(s+k*I)%*%s-d*solve(s+k*I)-I)%*%bb
bibets<-bibet%*%t(bibet)
mse<-dbd+bibets
mse1<-sum(diag(mse))
return(mse1)
}
arr[i*j]<-altem1(formula,k[j],d[i],aa,data,na.action)
falt1=file("alt1.data","a+")
cat(k[j],d[i],arr[i*j],"\n",file=falt1,append=TRUE)
close(falt1)
}
}
mat<-read.table("alt1.data")
unlink("alt1.data")
rmat<-matrix(mat[,3L],nrow=NROW(d),dimnames=list(c(paste0("d=",d)),c(paste0("k=",k))))
return(rmat)
}
pl1<-plotalt1(formula,k,d,aa,data,na.action)
plotpralt1<-function(formula,k,d,aa,data=NULL,na.action,...)
{
j<-0
i<-0
arr<-0
for (j in 1:nrow(k))
{
for (i in 1:nrow(d))
{
apre<-function(formula,k,d,aa,data,na.action,...)
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
ab<-matrix(aa,NCOL(x))
a<-t(ab)
pr<-0
 i<-1
 m<-1
for (i in 1:nrow(x))
{
               subsum<-0
bb<-c(x[i,]%*%t(x[i,]))
 z<-(solve(s+k*I-bb)-d*solve(s+k*I-bb)%*%solve(s-bb))%*%(t(x)%*%y-x[i,]*y[i])
            
for (m in 1:ncol(x))      
subsum<-subsum+(a[m]%*%x[i,m]%*%z[m])         
pr<-pr+(y[i]-subsum)^2
}
pre<-t(pr)
return(pre)
}
arr[j*i]<-apre(formula,k[j],d[i],aa,data,na.action)
fpralte1=file("praltev1.data","a+")
cat(k[j],d[i],arr[j*i],"\n",file=fpralte1,append=TRUE)
close(fpralte1)
}
}
pmat<-read.table("praltev1.data")
unlink("praltev1.data")
rpmat<-matrix(pmat[,3L],nrow=NROW(d),dimnames=list(c(paste0("d=",d)),c(paste0("k=",k))))
return(rpmat)
}
pralt1<-plotpralt1(formula,k,d,aa,data,na.action)
if (!press) prmse<-pl1 else prmse<-pralt1
if (nrow(k)>1L | nrow(d)>1L) val<-prmse else val<-npt
val
}
