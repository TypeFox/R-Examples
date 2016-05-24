trim.test<-function(Y,X,tr=.2){
n.i<-tapply(Y,X,length)
r<-nlevels(X)
Xn<-as.numeric(X)
h<-matrix(nrow=r,ncol=1)
d<-matrix(nrow=r,ncol=1)
for(i in 1:r){
h[i]<-length(trim.me(Y[Xn==i],trim=tr))
h}
for(i in 1:r){
d[i]<-((n.i[i]-1)*var(win(Y[Xn==i])))/(h[i]*(h[i]-1))
d
}
w<-1/d
U<-sum(w)
X.tilda<-(1/U)*sum(w*as.matrix(tapply(Y,X,function(x){mean(x,tr)})))
A<-(1/(r-1))*sum(w*as.matrix((tapply(Y,X,function(x){mean(x,tr)})-X.tilda)^2))
B<-(2*(r-2)/(r^2-1))*sum(((1-w/U)^2)/h)
F.star<-A/(1+B)
nu1<-r-1
nu2<-1/((3/(r^2-1))*sum(((1-(w/U))^2)/(h-1)))
res<-list()
res$Results<-data.frame(nu1=nu1,nu2=nu2,F=F.star,P.val=pf(F.star,nu1,nu2,lower.tail=FALSE))
colnames(res$Results)<-c("df1","df2","F*","P(>F)")
res
}