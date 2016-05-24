simul.lmp<-function(n, m, q, data, int=0, sigmap=1, p=2, lp=FALSE){
if (!is.numeric(n)||!is.numeric(m)||!is.numeric(q)||!is.numeric(data)||!is.numeric(int)
||!is.numeric(sigmap)||!is.numeric(p)) 
stop("Non-numeric argument to mathematical function")
if (length(data)!=q) stop("Length of data vector must be q")
name<-paste("x",1:q,sep="")
ris<-matrix(nrow=m,ncol=q+3,dimnames=list(c(1:m),c("(intercept)",name,"Sp","p")))
conv<-vector(length=m)
fr<-"y~"
for(k in 1:q) {fr<-paste(fr,name[k],sep="+")}
frm<-as.formula(fr)
if (lp==FALSE) pp<-NULL
if (lp==TRUE)  pp<-p
for (i in 1:m){
e<-rnormp(n,mu=0,sigmap=sigmap,p=p)
X<-matrix(runif(q*n),nrow=n,ncol=q,dimnames=(list(c(1:n),name)))
y<- int + data %*% t(X) + e
y<-as.vector(y)
X<-as.data.frame(X)
ll<-lmp(formula=frm,data=as.list(X),p=pp)
ll<-summary(ll)
ris[i,]<-c(ll$coeff,ll$sigma,ll$p)
if(lp==FALSE) conv[i]<-c(ll$iter)
}
mat<-ris
  ma<-matrix(nrow=2,ncol=ncol(ris),dimnames=list(c("Mean","Variance"),colnames(ris)))
  for(j in 1:ncol(ma))  {
  ma[1,j]<-mean(mat[,j])
  ma[2,j]<-(m-1)*var(mat[,j])/m
  }
ret<-list(dat=ris,table=ma)
ret$p<-p
ret$par<-c(int,data)
names(ret$par)<-colnames(ma)[1:length(ret$par)]
ret$frm<-frm
if(lp==FALSE) ret$iter<-sum(conv)
ret$lp<-lp
class(ret)<-c("simul.lmp","simul.mp")
ret
}

