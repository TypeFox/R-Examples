tail.plot<-function(dendat,type="both",split=median(dendat),
col="black",denmat=NULL,paletti=NULL,xlim=NULL,cex.axis=1,
pch=20,pchs=rep(20,1000),log="y")
{

if (is.null(paletti))
 paletti<-c("red","blue","green",
 "orange","navy","darkgreen",
 "orchid","aquamarine","turquoise",
 "pink","violet","magenta","chocolate","cyan",
 colors()[50:657],colors()[50:657])

if (type=="right.tail"){
   redu.ind<-(dendat>split) 
   dendat.redu<-dendat[redu.ind]
   ordi<-order(dendat.redu)
   dendat.ord<-dendat.redu[ordi]
   nredu<-length(dendat.redu)
   level<-seq(nredu,1)
   plot(dendat.ord,level,log=log,xlab="",ylab="",cex.axis=cex.axis,xlim=xlim,
   pch=pch)

   if (!is.null(denmat)){
      lkm<-dim(denmat)[2]
      for (i in 1:lkm){
          dencur<-denmat[,i]
          split=median(dencur)
          redu.ind<-(dencur>split) 
          dendat.redu<-dencur[redu.ind]
          ordi<-order(dendat.redu)
          dendat.ord<-dendat.redu[ordi]
          nredu<-length(dendat.redu)
          level<-seq(nredu,1)
          matplot(dendat.ord,level,log=log,xlab="",ylab="",cex.axis=cex.axis,
          add=TRUE,col=paletti[i],pch=pchs[i])
      }
   }
}

if (type=="left.tail"){
    redu.ind<-(dendat<split)
    dendat.redu<--dendat[redu.ind]
    ordi<-order(dendat.redu)
    dendat.ord<-dendat.redu[ordi]
    nredu<-length(dendat.redu)
    level<-seq(nredu,1)
    plot(-dendat.ord,level,log=log,xlab="",ylab="",cex.axis=cex.axis,xlim=xlim,
    pch=pch)

    if (!is.null(denmat)){
      lkm<-dim(denmat)[2]
      for (i in 1:lkm){
          dencur<-denmat[,i]
          split=median(dencur)
          redu.ind<-(dencur<split) 
          dendat.redu<--dencur[redu.ind]
          ordi<-order(dendat.redu)
          dendat.ord<-dendat.redu[ordi]
          nredu<-length(dendat.redu)
          level<-seq(nredu,1)
          matplot(-dendat.ord,level,log=log,xlab="",ylab="",cex.axis=cex.axis,
          add=TRUE,col=paletti[i],pch=pchs[i])
      }
    }
}

if (type=="both"){
    redu.ind<-(dendat<split)
    dendat.redu<--dendat[redu.ind]
    ordi<-order(dendat.redu)
    dendat.ord<-dendat.redu[ordi]
    nredu<-length(dendat.redu)
    level<-seq(nredu,1)

    redu.ind<-(dendat>split)
    dendat.redu<-dendat[redu.ind]
    ordi<-order(dendat.redu)
    dendat.ord.right<-dendat.redu[ordi]
    nredu<-length(dendat.redu)
    level.right<-seq(nredu,1)

    plot(c(-dendat.ord,dendat.ord.right),c(level,level.right),log=log,
    xlab="",ylab="",cex.axis=cex.axis,xlim=xlim,pch=pch)

    if (!is.null(denmat)){
      lkm<-dim(denmat)[2]
      for (i in 1:lkm){
          dencur<-denmat[,i]
          split=median(dencur)

          redu.ind<-(dencur<split) 
          dendat.redu<--dencur[redu.ind]
          ordi<-order(dendat.redu)
          dendat.ord<-dendat.redu[ordi]
          nredu<-length(dendat.redu)
          level<-seq(nredu,1)
          matplot(-dendat.ord,level,xlab="",ylab="",cex.axis=cex.axis,log=log,
          add=TRUE,col=paletti[i],pch=pchs[i])

          redu.ind<-(dencur>split) 
          dendat.redu<-dencur[redu.ind]
          ordi<-order(dendat.redu)
          dendat.ord<-dendat.redu[ordi]
          nredu<-length(dendat.redu)
          level<-seq(nredu,1)
          matplot(dendat.ord,level,xlab="",ylab="",cex.axis=cex.axis,log=log,
          add=TRUE,col=paletti[i],pch=pchs[i])
      }
    }
}

}











