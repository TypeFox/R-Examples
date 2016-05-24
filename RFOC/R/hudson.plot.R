hudson.plot<-function(m,col='red',pch=21,lwd=2,cex=1, bg='white' )
{
   m0=m2tk(m)
   uv=tk2uv(m0$T,m0$k)
   points(uv$u,uv$v,pch=pch,col=col,lwd=lwd,cex=cex, bg=bg)
}
