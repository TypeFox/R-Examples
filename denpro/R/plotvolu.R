plotvolu<-function(lst,length=NULL,
toplot=TRUE,data=FALSE,crit=NULL,orderrule="distcenter",
modelabel=FALSE,ptext=0,leimat=NULL,symbo=NULL,
info=NULL,infolift=0,infopos=0,
xmarginleft=0,xmarginright=0,ymargin=0,
xlim=NULL,ylim=NULL,
col="black",col.axis="black",
cutlev=NULL,xaxt="s",yaxt="s",
exmavisu=NULL,bg="transparent",tyyppi="n",
lty="solid",colo=FALSE,lowest="dens",proba=FALSE,
paletti=NULL,cex=NULL,modecolo=NULL,modepointer=NULL,upper=TRUE,
cex.axis=1,xlab="",ylab="",cex.lab=1,colothre=NULL,nodes=NULL)
{
if (upper) firstlevel<-min(lst$level) else firstlevel<-max(lst$level)
if (lowest=="dens") firstlevel<-0

parents<-lst$parent
levels<-lst$level
length<-lst$volume
if (proba) length<-lst$proba
center<-lst$center

mut<-multitree(parents)
if (is.null(lst$roots)) roots<-mut$roots else roots<-lst$roots
child<-mut$child
sibling<-mut$sibling

d<-dim(center)[1]
if (is.null(crit)){
   crit<-rep(0,d)          #order so that 1st closest to origo
   if (d==1) crit<-max(center)
   if (!is.null(lst$refe)) crit<-lst$refe
}

if (orderrule=="distcenter") sibord<-siborder(mut,crit,lst$distcenter)
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
  if (orderrule=="distcenter") sibord<-siborder(mut,crit,lst$distcenter)
  else sibord<-siborder(mut,crit,center)
  rootnum<-length(roots) 
  apuvecs<-matrix(NA,itemnum,4)
  for (i in 1:rootnum){
     inde<-roots[i]
     apuvecs[inde,]<-vecs[inde,]
     if (i==1) miniroot<-apuvecs[inde,1]
     else if (apuvecs[inde,1]<=miniroot) miniroot<-apuvecs[inde,1]
  }
  vecs<-apuvecs          #we give for the roots the previous positions
  vecs<-plotdata(roots,child,sibling,sibord,levels,length,vecs)
}

#####################################

depths<-NULL
segme<-T
lift<-NULL
modetest<-NULL
alpha<-NULL
axes<-T
modecolors<-NULL
modethickness<-1
leafcolors<-NULL
leaflift<-0
leafsymbo<-20
modelabels<-NULL
log<-""

nodenum<-length(vecs[,1])
xcoor<-matrix(0,2*nodenum,1)
ycoor<-matrix(0,2*nodenum,1)

for (i in 1:nodenum){
 xcoor[2*i-1]<-vecs[i,1]
 xcoor[2*i]<-vecs[i,3]
 ycoor[2*i-1]<-vecs[i,2]
 ycoor[2*i]<-vecs[i,4]
}

oriminnu<-min(orivecs[,1],na.rm=T)
minnu<-min(xcoor,na.rm=T)
if (is.null(cutlev)) xcoor<-xcoor-minnu
else xcoor<-xcoor-oriminnu

if (lowest=="dens") lowesti<-0 else lowesti<-min(lst$level)
#xlim<-c(min(vecs[,1],na.rm=T)-xmarginleft,max(vecs[,3],na.rm=T)+xmarginright)
if (is.null(ylim)){
    ylim<-c(lowesti,max(ycoor,na.rm=T)+ptext+ymargin)
    if (!is.null(cutlev)) 
    ylim<-c(cutlev,max(ycoor,na.rm=T)+ptext+ymargin)
}

if (toplot){
par(bg=bg)
plot(xcoor[order(xcoor)],ycoor[order(xcoor)],  #xcoor,ycoor,
xlab=xlab,ylab=ylab,axes=axes,xlim=xlim,ylim=ylim,xaxt=xaxt,
col=col,col.axis=col.axis,yaxt=yaxt,log=log,
type=tyyppi,lty=lty,cex.axis=cex.axis,cex.lab=cex.lab)
}
###########################################################

if ((tyyppi=="n") && (toplot)){

thick<-1
col<-col #"black"

for (i in 1:nodenum){
if (!is.na(ycoor[2*i-1])){

    yc<-ycoor[2*i-1]

    pare<-parents[i]
    if (pare==0) lowlev<-firstlevel else lowlev<-levels[pare]

    segments(xcoor[2*i-1],lowlev,xcoor[2*i-1],yc,col=col,lwd=thick)
    segments(xcoor[2*i],lowlev,xcoor[2*i],yc,col=col,lwd=thick)

    if (child[i]==0){  #we are in leaf

       segments(xcoor[2*i-1],yc,xcoor[2*i],yc,col=col,lwd=thick)

    }
    else{

       yc<-ycoor[2*i-1]

       childnum<-1
       curchi<-child[i]
       while (sibling[curchi]!=0){
           curchi<-sibling[curchi]
           childnum<-childnum+1
       }

       sibpointer<-matrix(0,childnum,1)
       curchi<-child[i]
       sibpointer[sibord[curchi]]<-curchi
       while (sibling[curchi]!=0){
           curchi<-sibling[curchi]
           sibpointer[sibord[curchi]]<-curchi
       }

       curchi<-sibpointer[1]
       x1<-xcoor[2*curchi-1]      
       segments(xcoor[2*i-1],yc,x1,yc,col=col,lwd=thick)
       x0<-xcoor[2*curchi] 

       cn<-2
       while (cn<=childnum){
             curchi<-sibpointer[cn]
             x1<-xcoor[2*curchi-1] 
             segments(x0,yc,x1,yc,col=col,lwd=thick)
             x0<-xcoor[2*curchi] 
             cn<-cn+1
       }

       segments(x0,yc,xcoor[2*i],yc,col=col,lwd=thick)

    }
}
}

for (i in 1:nodenum){
   if (is.null(cutlev)){
     orivecs[i,1]<-orivecs[i,1]-minnu
     orivecs[i,3]<-orivecs[i,3]-minnu
   }
   else{
     orivecs[i,1]<-orivecs[i,1]-oriminnu
     orivecs[i,3]<-orivecs[i,3]-oriminnu
   }
}   
if (modelabel) 
modelab<-plottext(parents,orivecs,ptext,leimat,symbo=symbo,cex=cex)  


}  #tyyppi = "n"


if (!is.null(lst$predictor.node)) 
segments(
xcoor[2*lst$predictor.node-1],
ycoor[2*lst$predictor.node-1],
xcoor[2*lst$predictor.node],
ycoor[2*lst$predictor.node])



############################################# exmavisu start

if (colo) exmavisu<-roots #1

if (!is.null(exmavisu)){

if (colo){
  if (is.null(paletti))
    paletti<-c("red","blue","green",
    "orange","navy","darkgreen",
    "orchid","aquamarine","turquoise",
    "pink","violet","magenta","chocolate","cyan",
    colors()[50:657],colors()[50:657])

  col<-colobary(lst$parent,paletti,modecolo=modecolo,modepointer=modepointer)

  if (!is.null(colothre))
  col<-colobary.merge(lst$parent,lst$level,colothre,paletti)
  if (!is.null(nodes))
  col<-colobary.nodes(lst$parent,nodes,paletti)

}
else col<-rep("blue",length(lst$parent))

for (i in 1:length(exmavisu)){

node<-exmavisu[i]

x1<-xcoor[2*node-1] 
x2<-xcoor[2*node]
lev<-levels[node]
if (parents[node]>0) lev0<-levels[parents[node]] else lev0<-firstlevel
polygon(c(x1,x2,x2,x1),c(lev0,lev0,lev,lev),col=col[node],lty="blank")

pino<-matrix(0,nodenum,1)
pino[1]<-child[node]
if (child[node]>0) pinoin<-1 else pinoin<-0

while (pinoin>0){
   node<-pino[pinoin]
   pinoin<-pinoin-1   

   x1<-xcoor[2*node-1] 
   x2<-xcoor[2*node]
   lev<-levels[node]
   if (parents[node]>0) lev0<-levels[parents[node]] else lev0<-firstlevel
   polygon(c(x1,x2,x2,x1),c(lev0,lev0,lev,lev),col=col[node],lty="blank")

   if (sibling[node]>0){
         pinoin<-pinoin+1
         pino[pinoin]<-sibling[node] 
   }

   while (child[node]>0){    #go to left and put right nodes to stack
         node<-child[node]

         x1<-xcoor[2*node-1] 
         x2<-xcoor[2*node]
         lev<-levels[node]
         if (parents[node]>0) lev0<-levels[parents[node]] else lev0<-firstlevel
         polygon(c(x1,x2,x2,x1),c(lev0,lev0,lev,lev),col=col[node],lty="blank")

         if (sibling[node]>0){
            pinoin<-pinoin+1
            pino[pinoin]<-sibling[node] 
         }
   }
}
}
}
####################### exmavisu end

if (data) return(list(xcoor=xcoor,ycoor=ycoor))


}





