checkm <-
function(formula,data,na.action,...)
{
cal<-match.call(expand.dots=FALSE)
mat<-match(c("formula","data","na.action"), names(cal))
cal<-cal[c(1L,mat)]
cal[[1L]]<-as.name("model.frame") 
cal<-eval(cal)       
md<- attr(cal, "terms")
dmat<-model.matrix(md,cal,contrasts)
if (NCOL(dmat)==1L) svec<-dmat else svec<-dmat[,-1L]
smat<-as.matrix(svec)
cor_mat<-cor(smat)
eimat<-eigen(cor_mat)
cn<-sqrt(abs(max(eimat$values)/min(eimat$values)))
 invcor<-solve(cor(smat))
avecor<-t(diag(invcor))
res<-list("*****Condition number*****"=cn,"*****VIF*****"=avecor) 
res
}
