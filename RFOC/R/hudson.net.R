hudson.net<-function(add=FALSE, POINTS=TRUE, TEXT=TRUE, colint='grey', colext='black')
{
   ###  plot a Hudson Net

  if(missing(colint)) colint='grey'
 if(missing(colext))   colext='black'

  
   p1=c(0,1)
   p2=c(-4/3,-1/3)
   p3=c(0,-1)
   p4=c(4/3,1/3)

   bmar=0.1
   b0=4/3+bmar
   bwid=2;lwid=1   

   if(add==FALSE)
     {
   ## par(omi=c(0,0,0,0),mai=c(0,0,0,0),font=2)
   plot(b0*c(-1,1),b0*c(-1,1),axes=FALSE,type='n',main='',xlab='',ylab='')
 }

   ks=seq(from=-0.9,to=0.9,by=0.1)
   Ts=seq(from=-0.9,to=0.9,by=0.1)
   for(i in 1:length(Ts))
   {
      ks0=seq(from=-1,to=1,by=0.001)
      tau0=0
      v0=0;u0=0
      for(j in 1:length(ks0))
      {
         uv0=tk2uv(Ts[i],ks0[j])
         v0[j]=uv0$v;u0[j]=uv0$u
      }
      lines(u0,v0,lty=2,lwd=lwid,col=colint)
   } 

   for(i in 1:length(ks))
   {
      Ts0=seq(from=-1,to=1,by=0.001)
      v0=0;u0=0
      for(j in 1:length(Ts0))
      {
         uv0=tk2uv(Ts0[j],ks[i])
         v0[j]=uv0$v;u0[j]=uv0$u
      }
      lines(u0,v0,lty=2,lwd=lwid,col=colint)
   } 

   hp1=tk2uv(1,0);hp2=tk2uv(-1,0)
   vp1=tk2uv(0,1);vp2=tk2uv(0,-1)

   lines(c(hp1$u,hp2$u),c(hp1$v,hp2$v),col=colint,lty=1,lwd=bwid)
   lines(c(vp1$u,vp2$u),c(vp1$v,vp2$v),col=colint,lty=1,lwd=bwid)

   lines(c(p1[1],p2[1]),c(p1[2],p2[2]),col=colext,lty=1,lwd=bwid)
   lines(c(p2[1],p3[1]),c(p2[2],p3[2]),col=colext,lty=1,lwd=bwid)
   lines(c(p3[1],p4[1]),c(p3[2],p4[2]),col=colext,lty=1,lwd=bwid)
   lines(c(p4[1],p1[1]),c(p4[2],p1[2]),col=colext,lty=1,lwd=bwid)

   Volp=tk2uv(0,1)
   DCp=tk2uv(0,0);CLVDp=tk2uv(1,0)
   LVDp=tk2uv(-1,1/3);CKp=tk2uv(-1,5/9)
    
if(POINTS)
  {
   points(LVDp$u,LVDp$v,pch=8)
   
   points(-LVDp$u,-LVDp$v,pch=8)
   
   points(CKp$u,CKp$v,pch=8)
  
   points(-CKp$u,-CKp$v,pch=8)
 }

   if(TEXT)
     {
       text(DCp$u,DCp$v,'DC',font=2)
       text(CLVDp$u,CLVDp$v,pos=4,offset=0.4,'-CLVD',fon=2)
       text(-CLVDp$u,CLVDp$v,pos=2,offset=0.4,'+CLVD',fon=2)
       text(Volp$u,-Volp$v,pos=1,offset=0.4,'-V',fon=2)
       text(Volp$u,Volp$v,pos=3,offset=0.4,'+V',fon=2)
       text(LVDp$u,LVDp$v,pos=2,offset=0.4,'+LVD',fon=2)
       text(-LVDp$u,-LVDp$v,pos=4,offset=0.4,'-LVD',fon=2)
       text(CKp$u,CKp$v,pos=2,offset=0.4,'+Crack',fon=2)
       text(-CKp$u,-CKp$v,pos=4,offset=0.4,'-Crack',fon=2)
     }

   
}
