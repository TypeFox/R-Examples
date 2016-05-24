rls <-
function(formula,r,R,delt,data,na.action,...)
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
colnames(rb)<-c("Estimate")
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)    
dbd<-ev*solve(s)-ev*solve(s)%*%t(RR)%*%solve(RR%*%solve(s)%*%t(RR))%*%RR%*%solve(s)
Standard_error<-sqrt(diag(abs(dbd)))
dbd<-ev*solve(s)-ev*solve(s)%*%t(RR)%*%solve(RR%*%solve(s)%*%t(RR))%*%RR%*%solve(s)
bibet<-solve(s)%*%t(RR)%*%solve(RR%*%solve(s)%*%t(RR))%*%del
bibets<-bibet%*%t(bibet)
mse<-dbd+bibets
rdel<-matrix(del,nrow(RR))
lenr<-length(RR)
dlpt<-diag(RR%*%xin%*%t(RR))
if (lenr==ncol(RR)) ilpt<-sqrt(solve(abs(dlpt))) else ilpt<-sqrt(solve(diag(abs(dlpt))))
upt<-RR%*%rb
tb<-t(upt)
t_statistic<-((tb-t(rdel))%*%ilpt)/sqrt(ev)
tst<-t(2L*pt(-abs(t_statistic),df<-(NROW(x)-NCOL(x))))
pvalue<-c(tst,rep(NA,(NCOL(x)-NROW(RR))))
mse1<-sum(diag(mse))
mse1<-round(mse1, digits <- 4L)
names(mse1)<-c("MSE")
t_statistic<-c(t_statistic,rep(NA,(NCOL(x)-NROW(RR))))
ans1<-cbind(rb,Standard_error,t_statistic,pvalue)
ans<-round(ans1, digits <- 4L)
anw<-list("*****Restricted Least Square Estimator*****"=ans,"*****Mean square error value******"=mse1)
anw
}
