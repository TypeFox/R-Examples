

musicNMR = function(source,destination){
   FID=paste(source,"fid",sep="")
   zz <- file(FID, "rb"); 
   acqu=paste(source,"acqu",sep="")
   ss=scan(acqu, what="character", sep=NULL)
   time=round(as.numeric(ss[1+which(ss=="##$SW_h=")])*2)
   a1=readBin(zz, integer(), 256*1024);
   savewav(a1, f=time, filename=destination);
   close(zz)
}


plotFID=function(x,ADD=FALSE, ...){
  if(ADD==FALSE){
  plot(x,xlab="Time (seconds)",ylab="Intensity (a.u.)",type="l", ...)}
 else{
 points(x,type="l", ...)
 }
}


musicMatrix= function(ma,destination){
   x=ma[,1]
   time=length(x)/x[length(x)]
   a1=ma[,2]
   savewav(a1, f=time, filename=destination);
}
