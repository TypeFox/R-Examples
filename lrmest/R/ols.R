ols <-
function(formula,data,na.action,...)
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
colnames(bb) <- c("Estimate")
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
dbb<-ev*xin
Standard_error<-sqrt(diag(abs(dbb)))
dbt<-t(bb)
sdbd_inv<-(sqrt(diag(abs(dbb))))^-1
sdbd_inv_mat<-diag(sdbd_inv)
if (NCOL(dbt)==1L) tbb<-dbt*sdbd_inv else tbb<-dbt%*%sdbd_inv_mat
tst<-t(tbb)
pval<-t(2*pt(-abs(tbb),df<-(NROW(x)-NCOL(x))))
        colnames(pval) <- c("p_value")
colnames(tst) <- c("t_statistic")
mse1<-sum(diag(dbb))
names(mse1)<-c("MSE")
mse1<-round(mse1, digits <- 4L)
if (NCOL(x)==1L) svec<-x else svec<-x[,-1L]
smat<-as.matrix(svec)
         invcor<-solve(cor(smat))
avecor<-sum(diag(invcor))/NCOL(x[,-1L])
con<-if (avecor>10L) warning ("There is a multicollinearity") else "There is not multicollinearity" 
ans<-cbind(bb,Standard_error,tst,pval)
ans1<-round(ans, digits <- 4L)
adw<-list("*****Ordinary Least Square Estimator******"=ans1,"*****Mean square error value*****"=mse1)
adw
}
