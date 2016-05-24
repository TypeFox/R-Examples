plottree<-function(lst,
plot=T,data=F,crit=NULL,orderrule="distcenter",
modelabel=TRUE,ptext=0,leimat=NULL,symbo=NULL,
info=NULL,infolift=0,infopos=0,infochar=NULL,
xmarginleft=0,xmarginright=0,ymargin=0,
xlim=NULL,ylim=NULL,
col="black",col.axis="black",linecol=rep("black",length(lst$parent)),
pch=21,dimen=NULL,yaxt="s",axes=T,
cex=NULL,nodemag=NULL,linemag=1,cex.axis=1,ylab="",cex.lab=1,
colo=FALSE,paletti=NULL,lowest="dens")
{ 
# create vector verticalPos
# find modes, number of modes, attach vertical position to modes
# position of parent is the mid of positions of children:
# use multitree to find siblings of node and "parent" to fine parent
#
#pch=19: solid circle, pch=20: bullet (smaller circle), 
#pch=21: circle, pch=22: square, 
#pch=23: diamond, pch=24: triangle point-up, 
#pch=25: triangle point down. 

if (colo){
  if (is.null(paletti))
    paletti<-c("red","blue","green",
    "orange","navy","darkgreen",
    "orchid","aquamarine","turquoise",
    "pink","violet","magenta","chocolate","cyan",
    colors()[50:657],colors()[50:657])

  col<-colobary(lst$parent,paletti,modecolo=NULL,modepointer=NULL)
  linecol<-col
}
else col<-rep(col,length(lst$parent))

parent<-lst$parent
level<-lst$level
center<-lst$center
if (is.null(center)){
   nodenum<-length(parent)
   dimen<-length(lst$refe)
   nodenum<-length(lst$parent)
   center<-matrix(1,dimen,nodenum)
}
#      
mut<-multitree(parent)    #create multitree 
roots<-mut$roots
child<-mut$child
sibling<-mut$sibling 

if (is.null(dimen)){
  d<-dim(center)[1]
}
else{
  d<-dimen
}

if (is.null(crit)){
   crit<-rep(0,d)          #order so that 1st closest to origo
   if (d==1) crit<-max(center)
   if (!is.null(lst$refe)) crit<-lst$refe
}
if (orderrule=="distcenter") sibord<-siborder(mut,crit,lst$distcenter)
else sibord<-siborder(mut,crit,center)

mlkm<-moodilkm(parent)
modloc<-mlkm$modloc   
#mlkm$modnodes
modenum<-mlkm$lkm  

lst$center<-center
modelinks<-siborToModor(lst)        #make links in right order

itemnum<-length(parent)    
verticalPos<-matrix(0,itemnum,1)

step<-1/modenum
curloc<-0
for (i in 1:modenum){
   curmode<-modelinks[i]   
   verticalPos[curmode]<-curloc
   curloc<-curloc+step
} 


for (i in 1:modenum){
   curnode<-modloc[i]
   par<-parent[curnode]
   while (par>0){
      #calculate mid of children of par
      #go to the end of sibling list
        chi<-child[par]
        summa<-verticalPos[chi]
        childNum<-1
        while(sibling[chi]>0){
           chi<-sibling[chi]
           summa<-summa+verticalPos[chi]
           childNum<-childNum+1
        }                            
        verticalPos[par]<-summa/childNum
        par<-parent[par]
   }
}

if (lowest=="dens") lowesti<-0 else lowesti<-min(lst$level)
if (is.null(ylim)) ylim<-c(lowesti-ymargin,max(level)+ptext+ymargin)
xlim<-c(min(verticalPos)-xmarginleft,max(verticalPos)+xmarginright)
#axes<-
plot(verticalPos,level,xlab="",ylab=ylab,xlim=xlim,ylim=ylim,xaxt="n",
col=col,col.axis=col.axis,pch=pch,yaxt=yaxt,axes=axes,cex=nodemag,
cex.axis=cex.axis,cex.lab=cex.lab)  

for (i in 1:itemnum){
    if (parent[i]>0){
        xchild<-verticalPos[i]
        ychild<-level[i]
        xparent<-verticalPos[parent[i]]
        yparent<-level[parent[i]]
        segments(xparent,yparent,xchild,ychild,col=linecol[i],lwd=linemag)
     }
}                
#
# lets plot info
#
if (!is.null(info)){
   nodenum<-itemnum
   infolocx<-matrix(nodenum,1)
   infolocy<-matrix(nodenum,1)
   #
   for (i in 1:nodenum){
     infolocx[i]<-verticalPos[i] 
     infolocy[i]<-level[i]+infolift
   }
   digits<-3
   info<-format(info,digits=digits)
   adj<-NULL
   pos<-infopos
   text(infolocx,infolocy,info,pos,adj,cex=cex)       
}
#
# lets plot character info
#
if (!is.null(infochar)){
   nodenum<-itemnum
   infolocx<-matrix(nodenum,1)
   infolocy<-matrix(nodenum,1)
   #
   for (i in 1:nodenum){
     infolocx[i]<-verticalPos[i] 
     infolocy[i]<-level[i]+infolift
   }
   pos<-infopos
   text(infolocx,infolocy,infochar,pos,cex=cex)       
}
#
# lets plot labels for modes
#
if (modelabel){
#
xcoor<-verticalPos
ycoor<-level
#
mlkm<-moodilkm(parent)
modloc<-mlkm$modloc  
modenum<-length(modloc)
modelocx<-matrix(0,modenum,1)
modelocy<-matrix(0,modenum,1)
if (is.null(leimat)){
   if (is.null(symbo)){
       labels<-paste("M",1:modenum,sep="")
   }
   else{
      labels<-paste(symbo,1:modenum,sep="")
   }
}
else{
   labels<-leimat
}
xcor<-matrix(0,modenum,1)                       
for (i in 1:modenum){
    loc<-modloc[i]
    xcor[i]<-xcoor[loc]
}
modloc<-omaord2(modloc,xcor)
for (i in 1:modenum){
    loc<-modloc[i]
    modelocx[i]<-xcoor[loc]
    modelocy[i]<-ycoor[loc]+ptext
}
text(modelocx,modelocy,labels,cex=cex)      
##
}
###############
}









