profkern<-function(dendat,h,N,Q,cvol=TRUE,ccen=TRUE,cfre=FALSE,kernel="epane",
compoinfo=FALSE,trunc=3,threshold=0.0000001,sorsa="crc",hw=NULL)
{

if (kernel=="gauss") h<-h*trunc   #trunc<-3

hnum<-length(h)
hrun<-1
while (hrun<=hnum){
   hcur<-h[hrun]

   if (sorsa=="crc")
   curtree<-profkernCRC(dendat,hcur,N,Q,kernel=kernel,compoinfo=compoinfo,
            trunc=trunc,threshold=threshold,hw=hw)
   else
   curtree<-profkernC(dendat,hcur,N,Q)

   if (hrun==1){
      if (hnum==1){
          treelist<-curtree
      }
      else{
          treelist=list(curtree)
      }
   }
   else{
      treelist=c(treelist,list(curtree))
   }
   hrun<-hrun+1
}
#
return(treelist)
}
