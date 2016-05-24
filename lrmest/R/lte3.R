lte3 <-
function(formula,k,d,press=FALSE,data=NULL,na.action,...)
{
k<-as.matrix(k)
d<-as.matrix(d)
k1<-k[1L]
d1<-d[1L]
ltes3<-function(formula,k1,d1,data=NULL,na.action,...)
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
bk<-(solve(s+I)+d1*solve(s+I)%*%solve(s+k1*I))%*%t(x)%*%y
colnames(bk) <- c("Estimate")
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
dbd<-ev*(solve(s+I)+d1*solve(s+I)%*%solve(s+k1*I))%*%s%*%(solve(s+I)+d1*solve(s+k1*I)%*%solve(s+I))
Standard_error<-sqrt(diag(abs(dbd)))
dbt<-t(bk)
dbd<-ev*(solve(s+I)+d1*solve(s+I)%*%solve(s+k1*I))%*%s%*%(solve(s+I)+d1*solve(s+k1*I)%*%solve(s+I))
sdbd_inv<-(sqrt(diag(abs(dbd))))^-1
sdbd_inv_mat<-diag(sdbd_inv)
if (NCOL(dbt)==1L) tbd<-dbt*sdbd_inv else tbd<-dbt%*%sdbd_inv_mat
hggh<-t(tbd)
tst<-t(2L*pt(-abs(tbd),df<-(NROW(x)-NCOL(x))))
        colnames(tst) <- c("p_value")
colnames(hggh) <- c("t_statistic")
bibet<-((solve(s+I)+d1*solve(s+I)%*%solve(s+k1*I))%*%s-I)%*%bb
bibets<-bibet%*%t(bibet)
mse<-dbd+bibets
mse1<-sum(diag(mse))
mse1<-round(mse1, digits = 4L)
names(mse1)<-c("MSE")
pr<-0
 i<-1
 m<-1
for (i in 1:nrow(x))
{
               subsum<-0
bb<-c(x[i,]%*%t(x[i,]))
 z<-(solve(s+I-bb)+d1*solve(s+I-bb)%*%solve(s+k1*I-bb))%*%(t(x)%*%y-x[i,]*y[i])
            
for (m in 1:ncol(x))      
subsum<-subsum+(x[i,m]%*%z[m])         
pr<-pr+(y[i]-subsum)^2
}
pre<-t(pr)
u1<-c(pre,rep(NA,NCOL(x)-1))
pres<-matrix(u1,NCOL(x))
colnames(pre)<-c("PRESS")
ans1<-cbind(bk,Standard_error,hggh,tst)
ans<-round(ans1, digits = 4L)
pre<-round(pre, digits = 4L)
anw<-list("*****Type (3) Liu Estimator*****"=ans,"*****Mean Square Error value*****"=mse1,"*****Prediction Sum of Squares value*****"=pre)
return(anw)
}
npt<-ltes3(formula,k1,d1,data,na.action)
plotlt3<-function(formula,k,d,data=NULL,na.action,...)
{
j<-1
i<-0
arr<-0
for (j in 1:nrow(k))
{
for (i in 1:nrow(d))
{
ltem3<-function(formula,k,d,data,na.action,...)
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
dbd<-ev*(solve(s+I)+d*solve(s+I)%*%solve(s+k*I))%*%s%*%(solve(s+I)+d*solve(s+k*I)%*%solve(s+I))
bibet<-((solve(s+I)+d*solve(s+I)%*%solve(s+k*I))%*%s-I)%*%bb
bibets<-bibet%*%t(bibet)
mse<-dbd+bibets
mse1<-sum(diag(mse))
return(mse1)
}
arr[j*i]<-ltem3(formula,k[j],d[i],data,na.action)
flte3<-file("lte3v.data","a+")
cat(k[j],d[i],arr[j*i],"\n",file=flte3,append=TRUE)
close(flte3)
}
}
mat<-read.table("lte3v.data")
unlink("lte3v.data")
rmat<-matrix(mat[,3L],nrow=NROW(d),dimnames=list(c(paste0("d=",d)),c(paste0("k=",k))))
return(rmat)
}
plte3<-plotlt3(formula,k,d,data,na.action)
plotprlt3<-function(formula,k,d,data=NULL,na.action,...)
{
j<-0
i<-0
arr<-0
for (j in 1:nrow(k))
{
for (i in 1:nrow(d))
{
pre11<-function(formula,k,d,data,na.action,...)
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
pr<-0
 i<-1
 m<-1
for (i in 1:nrow(x))
{
               subsum<-0
bb<-c(x[i,]%*%t(x[i,]))
 z<-(solve(s+I-bb)+d*solve(s+I-bb)%*%solve(s+k*I-bb))%*%(t(x)%*%y-x[i,]*y[i])
            
for (m in 1:ncol(x))      
subsum<-subsum+(x[i,m]%*%z[m])         
pr<-pr+(y[i]-subsum)^2
}
pre<-t(pr)
return(pre)
}
arr[j*i]<-pre11(formula,k[j],d[i],data,na.action)
fprlt3<-file("prlt3v.data","a+")
cat(k[j],d[i],arr[j*i],"\n",file=fprlt3,append=TRUE)
close(fprlt3)
}
}
pmat<-read.table("prlt3v.data")
unlink("prlt3v.data")
rprmat<-matrix(pmat[,3L],nrow=NROW(d),dimnames=list(c(paste0("d=",d)),c(paste0("k=",k))))
return(rprmat)
}
prlt3<-plotprlt3(formula,k,d,data,na.action)
if (!press) prmse<-plte3 else prmse<-prlt3
if (nrow(k)>1L | nrow(k)>1L) val<-prmse else val<-npt
val
}
