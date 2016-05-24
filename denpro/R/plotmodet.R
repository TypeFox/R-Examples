plotmodet<-function(mt,coordi=1,colot=NULL,
shift=0,xlim=NULL,xlab="",ylab="",
horilines=NULL,
symbo=20,loga=NULL,lty="dashed",
cex.axis=1,title=TRUE,cex.sub=1,cex.lab=1,
xaxt="s",yaxt="s")
{
epsi<-0.0000001
if (!is.null(horilines)) horilines<-mt$hseq[horilines]

if (is.null(loga)) 
   if (!is.null(mt$type)){
       if (mt$type=="greedy") loga<-"not"
       if (mt$type=="bagghisto") loga<-"not"
       if (mt$type=="carthisto") loga<-"not"
       if (mt$type=="kernel") loga<-"y"
   }
   else loga<-"y"

d<-dim(mt$xcoor)[2]

if (is.null(colot)){
    as<-mt$colot
}
else if (colot=="black"){
   lenni<-length(mt$ycoor)
   as<-matrix("black",lenni,1)
}
else as<-colot

if (d==1) xvec<-mt$xcoor
else xvec<-mt$xcoor[,coordi]
yvec<-mt$ycoor
len<-length(xvec)
for (i in 1:len){
  j<-i+1
  while (j<=len){
    if ((xvec[i]<=xvec[j]+epsi)&&(xvec[i]>=xvec[j]-epsi)&& 
        (yvec[i]<=yvec[j]+epsi)&&(yvec[i]>=yvec[j]-epsi)){
        #&&(as[i]!=as[j])){
             #xvec[j]<-xvec[j]+shift
             xvec[i]<-xvec[i]+shift
    }    
    j<-j+1
  }
}
if (loga=="y")
plot(xvec,yvec,col=as,xlim=xlim,xlab=xlab,ylab=ylab,pch=symbo,log=loga,
     cex.axis=cex.axis,cex.lab=cex.lab,xaxt=xaxt,yaxt=yaxt)   
else
plot(xvec,yvec,col=as,xlim=xlim,xlab=xlab,ylab=ylab,pch=symbo,
     cex.axis=cex.axis,cex.lab=cex.lab,xaxt=xaxt,yaxt=yaxt) 

if (title) title(sub=paste("coordinate",as.character(coordi)),cex.sub=cex.sub)


if (!is.null(horilines)){
  xmin<-min(xvec)
  xmax<-max(xvec)
  horilen<-length(horilines)
  for (i in 1:horilen){
    segments(xmin,horilines[i],xmax,horilines[i],lty=lty)
  }
}

itemnum<-length(mt$parent)
for (i in 1:itemnum){
    if (mt$parent[i]>0){
        xchild<-mt$xcoor[i,coordi]
        #if (loga=="y") ychild<-mt$ycoor[i] else 
        ychild<-mt$ycoor[i]
        xparent<-mt$xcoor[mt$parent[i],coordi]
        #if (loga=="y") yparent<-mt$ycoor[mt$parent[i]] else 
        yparent<-mt$ycoor[mt$parent[i]]
        collo<-mt$colot[i]  #mt$parent[i]]
        segments(xparent,yparent,xchild,ychild,col=collo)
     }
}

}

