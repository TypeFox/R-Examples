EMURA.Clayton <-
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
old=0.00001
repeat{
  score=sum(d)-sum( old/(r.grid-1+old)  )
  info=sum(  (r.grid-1)/(r.grid-1+old)^2  )
  new=old+score/info
  if(abs(old-new)<0.00001){break}
  old=new
}
alpha.cl=new

temp=(d.o==0)&(r.diag>=m^a)

##### Estimation of c, Fx and Sy #####
A=sum( (r.diag[temp]/m/sc.diag[temp])^(1-alpha.cl)-((r.diag[temp]-1)/m/sc.diag[temp])^(1-alpha.cl) )
prop.cl=( 1/(A+(1/m)^(1-alpha.cl)) )^(1/(1-alpha.cl))

A=(prop.cl*r.diag/m/sc.diag)^(1-alpha.cl)-(prop.cl*(r.diag-1)/m/sc.diag)^(1-alpha.cl)
A[r.diag<m^a]=0
Sy.cl=( 1-cumsum(A*dd.o) )^(1/(1-alpha.cl))
A[1]=(prop.cl/m)^(1-alpha.cl)
Fx.cl=( cumsum(A*(1-d.o)) )^(1/(1-alpha.cl))

Fx.cl=Fx.cl[d.o==0];Sy.cl=Sy.cl[d.o==1]

if(plotX==TRUE){
  plot(  sort(x.trunc),Fx.cl,type="s",xlim=c(min(x.trunc),max(x.trunc)),lwd=2,
       xlab="Time",ylab="Distribution function of X" )
}
if(plotY==TRUE){
  plot(  sort(z.trunc),Sy.cl,type="s",xlim=c(min(z.trunc),max(z.trunc)),lwd=2,
       xlab="Time",ylab="Survival function of Y"  )
}

list("alpha"=alpha.cl,"tau"=-(alpha.cl-1)/(alpha.cl+1),c=prop.cl,Fx=Fx.cl,Sy=Sy.cl)

}
