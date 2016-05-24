mixe <-
function(formula,r,R,dpn,delt,data,na.action,...)
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
RC<-matrix(R,NCOL(s))
RR<-t(RC)
if (is.matrix(R)) RR<-R else RR<-RR
if (length(dpn)==1L) shi<-dpn
else if (is.matrix(dpn)) shi<-dpn else shi<-diag(dpn)
de1<-as.matrix(delt)
bb<-xin%*%t(x)%*%y
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
w1<-solve(s/ev+t(RR)%*%solve(shi)%*%RR)
w2<-(t(x)%*%y)/ev+t(RR)%*%solve(shi)%*%r
bm<-w1%*%w2
colnames(bm) <- c("Estimate")
dbd<-w1
Standard_error<-sqrt(diag(abs(dbd)))
dbd<-w1
rdel<-matrix(delt,NROW(RR))
lenr<-length(RR)
dlpt<-diag(RR%*%xin%*%t(RR))
if (lenr==ncol(RR)) ilpt<-sqrt(solve(abs(dlpt))) else ilpt<-sqrt(solve(diag(abs(dlpt))))
upt<-RR%*%bm
tb<-t(upt)
t_statistic<-((tb-t(rdel))%*%ilpt)/sqrt(ev)
tst<-t(2L*pt(-abs(t_statistic),df=(NROW(x)-NCOL(x))))
pvalue<-c(tst,rep(NA,(NCOL(x)-NROW(RR))))
bibet<-xin%*%t(RR)%*%solve((shi/ev)+RR%*%xin%*%t(RR))%*%de1
bibets<-bibet%*%t(bibet)
mse<-dbd+bibets
mse1<-sum(diag(mse))
mse1<-round(mse1, digits <- 4L)
names(mse1)<-c("MSE")
t_statistic<-c(t_statistic,rep(NA,(NCOL(x)-NROW(RR))))
ans1<-cbind(bm,Standard_error,t_statistic,pvalue)
ans<-round(ans1, digits <- 4L)
anw<-list("*****Mixed  Regression Estimator*****"=ans,"*****Mean square error value*****"=mse1)
anw
}
