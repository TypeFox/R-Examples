plotbary<-function(lst,coordi=1,
plot=TRUE,data=FALSE,crit=NULL,orderrule="distcenter",
modelabel=FALSE,ptext=0,leimat=NULL,symbo=NULL,
info=NULL,infolift=0,infopos=0,
xmarginleft=0,xmarginright=0,ymargin=0,
xlim=NULL,ylim=NULL,xaxt="s",yaxt="s",
nodesymbo=20,col=NULL,col.axis="black",collines=NULL,paletti=NULL,
shift=0,shiftindex=NULL,
modlabret=FALSE,modecolo=NULL,modepointer=NULL,colometh="lst",
colothre=min(lst$level),lines=TRUE,wedge=FALSE,lty.wedge=2,
title=TRUE,titletext="coordinate",
cex=NULL,nodemag=NULL,cex.sub=1,cex.axis=1,newtitle=FALSE,cex.lab=1,
lowest="dens",subtitle=NULL)
{

parent<-lst$parent
center<-lst$center
level<-lst$level

if (is.null(paletti))
 paletti<-c("red","blue","green",
 "orange","navy","darkgreen",
 "orchid","aquamarine","turquoise",
 "pink","violet","magenta","chocolate","cyan",
 colors()[50:657],colors()[50:657])
if (is.null(col)) 
   if (colometh=="lst")
            col<-colobary(parent,paletti,
                 modecolo=modecolo,modepointer=modepointer)
   else col<-colobary.roots(lst$parent,lst$level,paletti=paletti,
                            colothre=colothre)

if (is.null(collines)) collines<-col

nodenum<-length(parent)
xcoordinate<-center[coordi,]

if (is.null(xlim))
xlim<-c(min(xcoordinate)-xmarginleft,max(xcoordinate)+xmarginright)
if (lowest=="dens") lowesti<-0 else lowesti<-min(lst$level)
ylim<-c(lowesti,max(level)+ptext+ymargin)

if (newtitle) xlab<-paste(titletext,as.character(coordi))
else xlab<-""
plot(xcoordinate,level,xlab=xlab,ylab="",
xlim=xlim,ylim=ylim,xaxt=xaxt,yaxt=yaxt,
pch=nodesymbo,col=col,col.axis=col.axis,cex=nodemag,
cex.axis=cex.axis,cex.lab=cex.lab) 
if (!is.null(subtitle)){ 
   title<-FALSE
   title(sub=subtitle,cex.sub=cex.sub)
}
if (title) title(sub=paste(titletext,as.character(coordi)),cex.sub=cex.sub)


if (lines){
   for (i in 1:nodenum){
       if (parent[i]>0){
           xchild<-xcoordinate[i]
           ychild<-level[i]
           xparent<-xcoordinate[parent[i]]
           yparent<-level[parent[i]]
           if (length(collines)>1) colli<-collines[i] else colli<-collines
           segments(xparent,yparent,xchild,ychild,col=colli) 
        }
   }
}

if (wedge){
  maxx<-max(xcoordinate)
  minx<-min(xcoordinate)
  righthigh<-maxx-lst$refe[coordi]
  lefthigh<-lst$refe[coordi]-minx
  segments(lst$refe[coordi],0,maxx,righthigh,lty=lty.wedge)
  segments(lst$refe[coordi],0,minx,lefthigh,lty=lty.wedge)
}

#########
#########
if (modlabret) modelabel<-TRUE
if (modelabel){

data<-plotprof(lst,plot=F,data=T,cutlev=NULL,ptext=NULL,info=NULL,
infolift=0,infopos=0,crit=crit,orderrule=orderrule)
vecs<-data$vecs
mlkm<-moodilkm(parent)
modloc<-mlkm$modloc 

nodenum<-length(vecs[,1])
xcoor<-matrix(0,2*nodenum,1)
ycoor<-matrix(0,2*nodenum,1)

for (i in 1:nodenum){
 xcoor[2*i-1]<-vecs[i,1]
 xcoor[2*i]<-vecs[i,3]
 ycoor[2*i-1]<-vecs[i,2]
 ycoor[2*i]<-vecs[i,4]
}                     
moodinum<-length(modloc)
modelocx<-matrix(0,moodinum,1)
modelocy<-matrix(0,moodinum,1)
if (is.null(leimat)){
   if (is.null(symbo)){
       labels<-paste("M",1:moodinum,sep="")
   }
   else{
         if (symbo=="empty") labels<-paste("",1:moodinum,sep="")
         else  labels<-paste(symbo,1:moodinum,sep="")
   }
}
else{
   labels<-leimat
}
xcor<-matrix(0,moodinum,1)
for (i in 1:moodinum){
    loc<-modloc[i]
    xcor[i]<-xcoor[2*loc-1]
}
modloc<-omaord2(modloc,xcor)
for (i in 1:moodinum){
    loc<-modloc[i]
    modelocx[i]<-xcoordinate[loc]
    modelocy[i]<-level[loc]+ptext
}
if (!is.null(shiftindex)) modelocx[shiftindex]<-modelocx[shiftindex]+shift
text(modelocx,modelocy,labels,cex=cex)         

if (modlabret){ 
   d<-dim(lst$center)[1]
   modelocat<-matrix(0,moodinum,d)
   for (i in 1:moodinum){
       loc<-modloc[i]
       modelocat[i,]<-lst$center[,loc]
   }   
   return(list(modelocat=modelocat,labels=labels))
}

}
############################################
}











