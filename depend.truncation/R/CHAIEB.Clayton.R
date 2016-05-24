CHAIEB.Clayton <-
function(x.trunc,z.trunc,d,a=1/10,plotX=TRUE,plotY=TRUE){

m=length(x.trunc)

#### Risk set on diagonal line ####
t=c(x.trunc,z.trunc)
t.o=t[order(t)]
d.o=c(rep(0,m),rep(1,m))[order(t)]
dd.o=c(rep(0,m),d)[order(t)]
r.diag=numeric(2*m)
for(i in 1:(2*m)){r.diag[i]=sum( (x.trunc<=t.o[i])&(z.trunc>=t.o[i]) )}

#### K-M for censoring ####
sc.diag=cumprod( 1-d.o*(1-dd.o)/r.diag*(r.diag>=m^a) )
sc=(sc.diag[d.o==1])[rank(z.trunc)]

#### Risk set on upper wedge ####
r.point=sckm.point=r.grid=sckm.grid=NULL
for(i in 1:m){
  x=x.trunc[i] ;z=z.trunc[i]
  if(d[i]==1){
    r.point=c(r.point,sum((x.trunc<=x)&(z.trunc>=z)))
    sckm.point=c(sckm.point,sc[i])
  }
  index=(x.trunc<=x)&(z.trunc>x)&(z.trunc<=z)&(d==1)
  num=sum(index)
  if(num>0){
    r.vec=numeric(num) ;sckm.vec<-numeric(num)
    for(j in 1:num){
      zz=z.trunc[index][j]
      r.vec[j]=sum( (x.trunc<=x)&(z.trunc>=zz) )
      sckm.vec[j]=(sc[index])[j]
    }
    r.grid=c(r.grid,r.vec) ;sckm.grid=c(sckm.grid,sckm.vec)
  }
}
pi.point=r.point/m ;pi.grid=r.grid/m

#### Estimation of alpha #####
x.orderx=x.trunc[order(x.trunc)] ;z.orderx=z.trunc[order(x.trunc)] ;d.orderx=d[order(x.trunc)]
con=dis=0
for(i in 2:m){
  x=x.orderx[i] ;z=z.orderx[i]
  for(j in 1:(i-1)){
    x12=x
    z12=min(z.orderx[j],z)
    if(x12<z12){
      dis=dis+d.orderx[i]*(z.orderx[j]>z) ;con=con+d.orderx[j]*(z.orderx[j]<z)
    }
  }
}
alpha.tau=dis/con

temp=(d.o==0)&(r.diag>=m^a)

##### Estimation of c, Fx and Sy #####
A=sum( (r.diag[temp]/m/sc.diag[temp])^(1-alpha.tau)-((r.diag[temp]-1)/m/sc.diag[temp])^(1-alpha.tau) )
prop.tau=( 1/(A+(1/m)^(1-alpha.tau)) )^(1/(1-alpha.tau))

A=(prop.tau*r.diag/m/sc.diag)^(1-alpha.tau)-(prop.tau*(r.diag-1)/m/sc.diag)^(1-alpha.tau)
A[r.diag<m^a]<-0
Sy.tau=( 1-cumsum(A*dd.o) )^(1/(1-alpha.tau))
A[1]=(prop.tau/m)^(1-alpha.tau)
Fx.tau=( cumsum(A*(1-d.o)) )^(1/(1-alpha.tau))

Fx.tau=Fx.tau[d.o==0];Sy.tau=Sy.tau[d.o==1]

if(plotX==TRUE){
  plot(  sort(x.trunc),Fx.tau,type="s",xlim=c(min(x.trunc),max(x.trunc)),lwd=2,
       xlab="Time",ylab="Distribution function of X" )
}
if(plotY==TRUE){
  plot(  sort(z.trunc),Sy.tau,type="s",xlim=c(min(z.trunc),max(z.trunc)),lwd=2,
       xlab="Time",ylab="Survival function of Y"  )
}

list("alpha"=alpha.tau,"tau"=-(alpha.tau-1)/(alpha.tau+1),c=prop.tau,Fx=Fx.tau,Sy=Sy.tau)

}
