# Chapter 5 R commands

# Section 5.2.1

data(eurodip)

year81=eurodip[,1]
nplus=sum(year81>0)
N=max(nplus,1)
rangd=N:(10^4*N)
post=1/(rangd*(rangd+1))
1/sum(post)
post=post/sum(post)
min(rangd[cumsum(post)>.5])

S=readline(prompt="Type  <Return>   to continue : ")

sum(pbino(nplus)*1:400)

# Section 5.2.2

S=readline(prompt="Type  <Return>   to continue : ")

prob=lchoose((471570:10^7),471570)+lgamma(2*(471570:10^7)-582681+1)-lgamma(2*(471570:10^7)+2)
range(prob)

S=readline(prompt="Type  <Return>   to continue : ")

n1=sum(eurodip[,1]>0)
n2=sum(eurodip[,2]>0)
m2=sum((eurodip[,1]>0) & (eurodip[,2]>0))
nc=n1+n2
nplus=n1+n2-m2

sum((1:400)*pcapture(2,nplus,nc))

S=readline(prompt="Type  <Return>   to continue : ")

sum((1:400)*pdarroch(n1,n2,m2))

S=readline(prompt="Type  <Return>   to continue : ")

for (i in 6:16) print(round(sum(pdarroch(n1,n2,i)*1:400)))

S=readline(prompt="Type  <Return>   to continue : ")

# Section 5.2.2

n3=sum(eurodip[,3]>0)
n4=sum(eurodip[,4]>0)
n5=sum(eurodip[,5]>0)
n6=sum(eurodip[,6]>0)
n7=sum(eurodip[,7]>0)
m3=sum(eurodip[,3]>0 & (eurodip[,1]>0 | eurodip[,2]>0))
m4=sum(eurodip[,4]>0 & (eurodip[,1]>0 | eurodip[,2]>0 | eurodip[,3]>0))
m5=sum(eurodip[,5]>0 & (eurodip[,1]>0 | eurodip[,2]>0 | eurodip[,3]>0 | eurodip[,4]>0))
m6=sum(eurodip[,6]>0 & (eurodip[,1]>0 | eurodip[,2]>0 | eurodip[,3]>0 | eurodip[,4]>0 | eurodip[,5]>0))
m7=sum(eurodip[,7]>0 & (eurodip[,1]>0 | eurodip[,2]>0 | eurodip[,3]>0 | eurodip[,4]>0 | eurodip[,5]>0 | eurodip[,6]))
nc=n1+n2+n3+n4+n5+n6+n7
nplus=nc-m2-m3-m4-m5-m6-m7
sum((1:400)*pcapture(7,nplus,nc))

S=readline(prompt="Type  <Return>   to continue : ")

S=2500;T=7;nplus=294;nc=519
lprob=lchoose(max(nplus,1):S,nplus)+
        lgamma(T*max(nplus,1):S-nc+1)-lgamma(T*max(nplus,1):S+2)
prob=c(rep(0,max(nplus,1)-1),exp(lprob-max(lprob)))
sum((1:S)*prob)/sum(prob)

S=readline(prompt="Type  <Return>   to continue : ")

lambda=200
nsimu=10^4
p=rep(1,nsimu); N=p
N[1]=2*nplus
p[1]=rbeta(1,nc+1,T*lambda-nc+1)
for (i in 2:nsimu){
  N[i]=nplus+rpois(1,lambda*(1-p[i-1])^T)
  p[i]=rbeta(1,nc+1,T*N[i]-nc+1)
  }

mean(N)
mean(p)
1/(1+2*sum(acf(N)$acf[-1]))
1/(1+2*sum(acf(p)$acf[-1]))

S=readline(prompt="Type  <Return>   to continue : ")

# Section 5.3

c2=m2
c3=sum(eurodip[,3]>0 & eurodip[,1]>0 )
g2=gibbscap1(10000,n1,c2,c3,200,10,5)

mean(g2$N)
mean(g2$p)
sq=seq(1,10000,by=20)
par(mfrow=c(5,2),mar=c(2,4,1,1))
plot(sq,g2$p[sq],ylab="p",xlab="",type="l")
hist(g2$p,nclass=100,prob=T,main="",ylab="")
plot(sq,g2$q[sq],ylab="q",xlab="",type="l")
hist(g2$q,nclass=100,prob=T,main="",ylab="")
plot(sq,g2$N[sq],ylab="N",xlab="",type="l")
hist(g2$N,nclass=100,prob=T,main="",ylab="")
plot(sq,g2$r1[sq],ylab=expression(r[1]),xlab="",type="l")
hist(g2$r1,nclass=100,prob=T,main="",ylab="")
plot(sq,g2$r2[sq],ylab=expression(r[2]),xlab="",type="l")
hist(g2$r2,nclass=100,prob=T,main="",ylab="")

S=readline(prompt="Type  <Return>   to continue : ")

dev.off()
plot(jitter(g2$r1,factor=2),jitter(g2$r2,factor=2),cex=0.5,
xlab=expression(r[1]),ylab=expression(r[2]),xlim=c(-0.5,5),ylim=c(-0.5,6),col=heat.colors(200))

S=readline(prompt="Type  <Return>   to continue : ")

# Section 5.4

mean(thresh(0:11,n1,c2,c3,1,0.1))
ar1=ardipper(10000,n1,c2,c3,1,0.1)
ar1=factor(ar1)
barplot(summary(ar1)/10000,main="",xlab=expression(r[1]))

S=readline(prompt="Type  <Return>   to continue : ")

# Section 5.5

solbeta(.1,.5,.3,10^(-4))

g3=gibbscap2(10000,eurodip)
apply(g3$p,2,mean)
apply(g3$phi,2,mean)
mean(g3$psi[3,3,])
mean(g3$psi[1,2,])

S=readline(prompt="Type  <Return>   to continue : ")

sq=seq(20,10000,by=20)
par(mfrow=c(3,2),mar=c(2,4,1,1))
plot(sq,g3$p[sq,1],ylab="p(1)",xlab="",type="l")
hist(g3$p[20:10000,1],nclass=100,prob=T,main="",ylab="")
plot(sq,g3$phi[sq,2],ylab=expression(phi(2)),xlab="",type="l")
hist(g3$phi[20:10000,2],nclass=100,prob=T,main="",ylab="")
plot(sq,g3$psi[3,3,sq],ylab=expression(psi(3,3)),xlab="",type="l")
hist(g3$psi[3,3,20:10000],nclass=100,prob=T,main="",ylab="")

