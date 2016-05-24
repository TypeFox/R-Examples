rota.seq<-function(dendat,col,pcf,ste=3,cut=dim(dendat)[1],simple=TRUE)
{
i1<-1
i2<-2
i3<-3
aa<-seq(0,2*pi,ste)
bb<-seq(0,2*pi,ste)
cc<-seq(0,2*pi,ste)
ii<-1
while (ii<=length(aa)){
  alpha<-aa[ii]
  jj<-1
  while (jj<=length(bb)){
    beta<-bb[jj]
    kk<-1
    while (kk<=length(cc)){
      gamma<-cc[kk]
      dexdat<-rotation3d(dendat,alpha,beta,gamma)
      plot.histdata(dexdat,col,pcf,i1,i2,i3,simple=simple,cut=cut,
      xlab=paste(as.character(ii),as.character(jj),as.character(kk)))
      kk<-kk+1
    }
  jj<-jj+1
  }
ii<-ii+1
}


}
