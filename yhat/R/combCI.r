combCI<-function(lowerCI,upperCI,est,star=FALSE){
  lowerCI<-format(lowerCI,3,format="f")
  upperCI<-format(upperCI,3,format="f")
  est<-format(est,3,format="f")
  combCI<-est
  for (i in 1:nrow(combCI)){
    for (j in 1: ncol(combCI)){
      if (combCI[i,j]!="NA"){
        combCI[i,j]<-paste(combCI[i,j],lowerCI[i,j],sep="(")
        combCI[i,j]<-paste(combCI[i,j],upperCI[i,j],sep=",")
        combCI[i,j]<-paste(combCI[i,j],")",sep="")
        if (star){
          if ((lowerCI[i,j] < 0) & (upperCI[i,j] > 0))
            combCI[i,j]<-paste(combCI[i,j]," ",sep="")
          else
            combCI[i,j]<-paste(combCI[i,j],"*",sep="")
        }
        combCI[i,j]<-toString(combCI[i,j])
       }
     }
   }
combCI<-as.data.frame(combCI)
return(combCI)
}