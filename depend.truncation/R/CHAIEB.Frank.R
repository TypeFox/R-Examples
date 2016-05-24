CHAIEB.Frank <-
function(x.trunc,z.trunc,d,a=1/10,plotX=TRUE,plotY=TRUE){

m=length(x.trunc)

#### Risk set on diagonal line ####
t=c(x.trunc,z.trunc)
t.o=t[order(t)]
d.o=c(rep(0,m),rep(1,m))[order(t)]
dd.o=c(rep(0,m),d)[order(t)]
r.diag=numeric(2*m)
for(i in 1:(2*m)){
r.diag[i]=sum( (x.trunc<=t.o[i])&(z.trunc>=t.o[i]) )
}

#### K-M for censoring ####
sc.diag=cumprod( 1-d.o*(1-dd.o)/r.diag*(r.diag>=m^a) )
sc=(sc.diag[d.o==1])[rank(z.trunc)]

#### Risk set on upper wedge ####
r.point=sckm.point=r.grid=sckm.grid=NULL
for(i in 1:m){
x=x.trunc[i]; z=z.trunc[i]
if(d[i]==1){
r.point=c(r.point,sum((x.trunc<=x)&(z.trunc>=z)))
sckm.point=c(sckm.point,sc[i])
}
index=(x.trunc<=x)&(z.trunc>x)&(z.trunc<=z)&(d==1)
num=sum(index)
if(num>0){
  r.vec=numeric(num); sckm.vec<-numeric(num)
  for(j in 1:num){
  zz=z.trunc[index][j]
  r.vec[j]=sum( (x.trunc<=x)&(z.trunc>=zz) ); sckm.vec[j]=(sc[index])[j]
  }
  r.grid=c(r.grid,r.vec); sckm.grid<-c(sckm.grid,sckm.vec)
}
}
pi.point=r.point/m; pi.grid=r.grid/m

#### Estimation of association #####
l=log(0.00001); u=log(10000)
repeat{
  mid=(l+u)/2
  theta.point=pi.point/sckm.point*l/( exp(pi.point/sckm.point*l)-1 )
  theta.grid=pi.grid/sckm.grid*l/( exp(pi.grid/sckm.grid*l)-1 )
  p0.point=r.point-1+theta.point
  score.l=sum( p0.point/(1+theta.point) )-sum( theta.grid/(1+theta.grid) )
  theta.point=pi.point/sckm.point*mid/( exp(pi.point/sckm.point*mid)-1 )
  theta.grid=pi.grid/sckm.grid*mid/( exp(pi.grid/sckm.grid*mid)-1 )
  p0.point=r.point-1+theta.point
  score.m=sum( p0.point/(1+theta.point) )-sum( theta.grid/(1+theta.grid) )
  if( score.l*score.m>0 ){l=mid}
  else{u=mid}
  if( abs(l-u)<0.0000001 ){break}
}
gamma.tau=(l+u)/2

##### Estimation of c, Fx and Sy #####
a.tau=a
repeat{
temp=(d.o==0)&(r.diag>=m^a.tau)
A=prod( (exp(gamma.tau*r.diag[temp]/m/sc.diag[temp])-1)/
( exp(gamma.tau*(r.diag[temp]-1)/m/sc.diag[temp])-1) )
alpha.tau=1+( exp(gamma.tau/m)-1 )*A
if(alpha.tau>0){break}
a.tau=a.tau*1.5 # increase the tuning parameter
}
prop.tau=gamma.tau/log(alpha.tau)

A=(alpha.tau^(prop.tau*r.diag/m/sc.diag)-1)/(alpha.tau^(prop.tau*(r.diag-1)/m/sc.diag)-1)
A[r.diag<m^a.tau]=1
Sy.tau=log( 1+(alpha.tau-1)/cumprod(A^dd.o),base=alpha.tau )
A=A[-1]
Fx.tau=c(prop.tau/m,log( 1+(alpha.tau^(prop.tau/m)-1)/cumprod( (1/A)^(1-d.o[-1]) ),base=alpha.tau ))

Fx.tau=Fx.tau[d.o==0]; Sy.tau=Sy.tau[d.o==1]

gamma=-log(alpha.tau)
debye=function(x){ x/(exp(x)-1)}
tau=-(1+4/gamma*(1/gamma*integrate(debye,0,gamma)$value-1))

if(plotX==TRUE){
  plot(  sort(x.trunc),Fx.tau,type="s",xlim=c(min(x.trunc),max(x.trunc)),lwd=2,
       xlab="Time",ylab="Distribution function of X" )
}
if(plotY==TRUE){
  plot(  sort(z.trunc),Sy.tau,type="s",xlim=c(min(z.trunc),max(z.trunc)),lwd=2,
       xlab="Time",ylab="Survival function of Y"  )
}

list("alpha"=1/alpha.tau,"tau"=tau,c=prop.tau,Fx=Fx.tau,Sy=Sy.tau)
#### 0<alpha<1(positive association); alpha>1(negative association) ####

}
