trim.ranef.test<-function(Y,X,tr=.2){
a<-nlevels(as.factor(X))
Xn<-as.numeric(X)

X.bar.Wi<-matrix(ncol=1,nrow=a)
S.sq.W.i<-matrix(ncol=1,nrow=a)
for(i in 1:a){
Y.ij<-win(Y[Xn==i])
X.bar.Wi[i]<-mean(Y.ij)
S.sq.W.i[i]<-var(Y.ij)
}

h<-matrix(nrow=a,ncol=1)
for(i in 1:a){
h[i]<-length(trim.me(Y[Xn==i],trim=tr))
h}

X.bar.ti<-tapply(Y,X,function(x){mean(x,trim=tr)})
X.bar.t<-mean(X.bar.ti)
BSST<-(1/(a-1))*sum((as.matrix(X.bar.ti)-X.bar.t)^2)

inner.WSSW<-matrix(nrow=a,ncol=1)
for(i in 1:a){
inner.WSSW[i]<-sum((win(Y[Xn==i]-X.bar.Wi[i])^2)/(h[i]*(h[i]-1)))
}

WSSW<-(1/a)*sum(inner.WSSW)
F.star<-BSST/WSSW

n.i<-tapply(Y,X,length)

d<-matrix(nrow=a,ncol=1)
for(i in 1:a){
d[i]<-((n.i[i]-1)*S.sq.W.i[i])/(h[i]*(h[i]-1))

}
nu1<-(((a-1)*sum(d))^2)/((sum(d))^2+(a-2)*a*sum(d^2))
nu2<-(sum(d))^2/(sum((d^2)/(h-1)))
res<-list()
res$Table<-data.frame(nu1=nu1,nu2=nu2,F=F.star,P.val=pf(F.star,nu1,nu2,lower.tail=FALSE))
colnames(res$Table)<-c("df1","df2","F*","P(>F)")
res
}