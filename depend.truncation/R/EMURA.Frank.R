EMURA.Frank <-
function(x.trunc,z.trunc,d,a=1/10,plotX=TRUE,plotY=TRUE){

m=length(x.trunc)

#### Risk set on diagonal line ####
t=c(x.trunc,z.trunc)
t.o=t[order(t)]
d.o=c(rep(0,m),rep(1,m))[order(t)]
dd.o=c(rep(0,m),d)[order(t)]
r.diag=numeric(2*m)
for(i in 1:(2*m)){r.diag[i]=sum( (x.trunc<=t.o[i])&(z.trunc>=t.o[i]))
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
  zz<-z.trunc[index][j]
  r.vec[j]=sum( (x.trunc<=x)&(z.trunc>=zz) ); sckm.vec[j]=(sc[index])[j]
  }
  r.grid=c(r.grid,r.vec); sckm.grid=c(sckm.grid,sckm.vec)
}
}
pi.point=r.point/m; pi.grid=r.grid/m

#### Estimation of association #####
old=log(0.00001)
repeat{
  theta.grid<-pi.grid/sckm.grid*old/(exp(pi.grid/sckm.grid*old)-1)
  p0.grid=r.grid-1+theta.grid
  w.grid=1-exp(pi.grid/sckm.grid*old)*theta.grid
  score=sum(d)-sum( theta.grid/p0.grid )
  info=sum( theta.grid*w.grid*(r.grid-1)/p0.grid^2 )/old
  new=old+score/info
  if(abs(new-old)<0.00001){break}
  old=new
}
gamma.lr=new

#### Quasi-Fisher scoring ####
old=gamma.lr
repeat{
  theta.point=pi.point/sckm.point*old/(exp(pi.point/sckm.point*old)-1)
  theta.grid=pi.grid/sckm.grid*old/(exp(pi.grid/sckm.grid*old)-1)
  p0.grid=r.grid-1+theta.grid
  w.point=1-exp(pi.point/sckm.point*old)*theta.point
  w.grid=1-exp(pi.grid/sckm.grid*old)*theta.grid
  score=sum(w.point)-sum(w.grid*theta.grid/p0.grid)
  info=sum( theta.grid*w.grid^2*(r.grid-1)/p0.grid^2 )/old
  new=old+score/info
  if(abs(new-old)<0.00001){break}
  old=new
}
gamma.cl=new

##### Estimation of c, Fx and Sy #####
a.cl=a
repeat{
  temp=(d.o==0)&(r.diag>=m^a.cl)
  A=prod( (exp(gamma.cl*r.diag[temp]/m/sc.diag[temp])-1)/
  ( exp(gamma.cl*(r.diag[temp]-1)/m/sc.diag[temp])-1) )
  alpha.cl=1+( exp(gamma.cl/m)-1 )*A
  if(alpha.cl>0){break}
  a.cl=a.cl*1.5  # increase the tuning parameter
}
prop.cl=gamma.cl/log(alpha.cl)

A=(alpha.cl^(prop.cl*r.diag/m/sc.diag)-1)/(alpha.cl^(prop.cl*(r.diag-1)/m/sc.diag)-1)
A[r.diag<m^a.cl]=1
Sy.cl=log( 1+(alpha.cl-1)/cumprod(A^dd.o),base=alpha.cl )
A=A[-1]
Fx.cl=c(prop.cl/m,log( 1+(alpha.cl^(prop.cl/m)-1)/cumprod( (1/A)^(1-d.o[-1]) ),base=alpha.cl ))

Fx.cl=Fx.cl[d.o==0]; Sy.cl=Sy.cl[d.o==1]

gamma=-log(alpha.cl)
debye=function(x){ x/(exp(x)-1)}
tau=-(1+4/gamma*(1/gamma*integrate(debye,0,gamma)$value-1))

if(plotX==TRUE){
  plot(  sort(x.trunc),Fx.cl,type="s",xlim=c(min(x.trunc),max(x.trunc)),lwd=2,
       xlab="Time",ylab="Distribution function of X" )
}
if(plotY==TRUE){
  plot(  sort(z.trunc),Sy.cl,type="s",xlim=c(min(z.trunc),max(z.trunc)),lwd=2,
       xlab="Time",ylab="Survival function of Y"  )
}

list("alpha"=1/alpha.cl,"tau"=tau,c=prop.cl,Fx=Fx.cl,Sy=Sy.cl)
#### 0<alpha<1(positive association); alpha>1(negative association) ####
  
}
