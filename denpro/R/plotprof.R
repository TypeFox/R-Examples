plotprof<-function(profile,length=NULL,
plot=TRUE,data=FALSE,crit=NULL,orderrule="distcenter",
modelabel=TRUE,ptext=0,leimat=NULL,symbo=NULL,
info=NULL,infolift=0,infopos=0,
xmarginleft=0,xmarginright=0,ymargin=0,
xlim=NULL,ylim=NULL,axes=TRUE,
col="black",col.axis="black",
cutlev=NULL,xaxt="n",exmavisu=NULL,cex.axis=1,cex=1)
{

#xaxs="e"    (extended)  not implemented?  xaxt="n"

parents<-profile$parent
levels<-profile$level
if (is.null(length)) length<-profile$volume
center<-profile$center

mut<-multitree(parents)
if (is.null(profile$roots)) roots<-mut$roots else roots<-profile$roots
child<-mut$child
sibling<-mut$sibling

d<-dim(center)[1]
if (is.null(crit)){
   crit<-rep(0,d)          #order so that 1st closest to origo
   if (d==1) crit<-max(center)
   if (!is.null(profile$refe)) crit<-profile$refe
}

if (orderrule=="distcenter") sibord<-siborder(mut,crit,profile$distcenter)
else sibord<-siborder(mut,crit,center)

itemnum<-length(parents)
vecs<-matrix(NA,itemnum,4)
vecs<-alloroot(vecs,roots,sibord,levels,length)
vecs<-plotdata(roots,child,sibling,sibord,levels,length,vecs)
orivecs<-vecs

if (!(is.null(cutlev))){
  cm<-cutmut(mut,cutlev,levels)              # cut the tree
  roots<-cm$roots
  sibling<-cm$sibling
  mut$roots<-roots
  if (orderrule=="distcenter") sibord<-siborder(mut,crit,profile$distcenter)
  else sibord<-siborder(mut,crit,center)
  rootnum<-length(roots) 
  apuvecs<-matrix(NA,itemnum,4)
  for (i in 1:rootnum){
     inde<-roots[i]
     apuvecs[inde,]<-vecs[inde,]
  }
  vecs<-apuvecs          #we give for the roots the previous positions
  vecs<-plotdata(roots,child,sibling,sibord,levels,length,vecs)
}

if (plot==TRUE){
   if (!(is.null(cutlev))){
     xlim<-c(omamin(vecs[,1])-xmarginleft,omamax(vecs[,3])+xmarginright)
     ylim<-c(omamin(vecs[,2]),omamax(vecs[,2])+ptext+ymargin)
   }
   else{
     xlim<-c(omamin(vecs[,1])-xmarginleft,omamax(vecs[,3])+xmarginright)
     if (is.null(ylim)) ylim<-c(0,omamax(vecs[,2])+ptext+ymargin)
   }
   plotvecs(vecs,segme=T,xlim=xlim,ylim=ylim,xaxt=xaxt,
   col=col,col.axis=col.axis,cex.axis=cex.axis,axes=axes)
   # use original vectors (numbering will be correct)
   if (modelabel){
      plottext(parents,orivecs,ptext,leimat,symbo,cex=cex)  
   }
   if (!is.null(info)){
      plotinfo(vecs,info,pos=infopos,adj=NULL,lift=infolift,digits=3)
   }
}
#
#
if (data==TRUE){
 return(list(sibord=t(sibord),vecs=vecs,parents=parents,levels=levels,
 length=length,center=center,remain=NULL))
}

}












