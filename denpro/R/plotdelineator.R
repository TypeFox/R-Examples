plotdelineator<-function(shtseq,coordi=1,ngrid=40,shift=0.05,
volumefunction=NULL,redu=TRUE,type="l")
{
if (is.null(volumefunction)){
   lnum<-length(shtseq$level)
   st<-shtseq$shtseq[[1]]
   td<-treedisc(st,shtseq$pcf,ngrid=ngrid)
   #td<-prunemodes(td,exmalim=0.5)$lst
   reduseq<-list(td)
   for (i in 2:lnum){
       st<-shtseq$shtseq[[i]]
       td<-treedisc(st,shtseq$pcf,ngrid=ngrid)
       #td<-prunemodes(td,exmalim=0.00001)$lst
       reduseq<-c(reduseq,list(td))
   }
   estiseq<-list(lstseq=reduseq,hseq=shtseq$level)
   mg<-modegraph(estiseq)
   plotmodet(mg,coordi=coordi,shift=shift)
}
else{
    vd<-volumefunction
    if (redu){
       x<-vd$delineator.redu[,coordi]
       y<-vd$delineatorlevel.redu
       or<-order(x)
       x1<-x[or]
       y1<-y[or]
       plot(x1,y1,type=type,
            ylab="level",xlab=paste("coordinate",as.character(coordi)))
    }
    else
       plot(vd$delineator[,coordi],vd$delineatorlevel,ylab="level")
}    

}


