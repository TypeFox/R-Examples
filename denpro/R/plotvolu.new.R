plotvolu.new<-function(lst,dens=TRUE)
{
mt<-multitree(lst$parent)
itemnum<-length(lst$volume)
rootnum<-length(mt$roots)
left<-mt$child
right<-mt$sibling
vecs<-matrix(0,itemnum,3)

sibord<-mt$siborder  #siborder.new(mt)

# allocate space for roots

rootsvolume<-0
for (i in 1:rootnum){
  now<-mt$roots[i]
  rootsvolume<-rootsvolume+lst$volume[now]
}
basis<-rootsvolume+rootsvolume/4
gaplen<-(basis-rootsvolume)/(rootnum+1)
rootlinks<-matrix(0,rootnum,1)  #make links in right order
{
if (rootnum==1){ 
  rootlinks[1]<-mt$roots[1]  #1
}
else{ 
     for (i in 1:rootnum){
         now<-mt$roots[i]
         roor<-sibord[now]
         rootlinks[roor]<-now
     }
}
xbeg<-0
xend<-0
for (i in 1:rootnum){
  now<-rootlinks[i]
  ycoo<-lst$level[now]
  xend<-xbeg+lst$volume[now]
  vecs[now,]<-c(xbeg,xend,ycoo)
  xbeg<-gaplen+xend
}
}
# allocate space for nonroots

for (i in 1:rootnum){
    pino<-matrix(0,itemnum,1)
    pino[1]<-mt$roots[i]  
    pinin<-1
    while (pinin>0){
        cur<-pino[pinin]      #take from stack
        pinin<-pinin-1
        if (left[cur]>0){     #if not leaf (root may be leaf)
           vecs<-allokoi.new(cur,vecs,lst,left,right,sibord)   
        }
        if (right[cur]>0){    #if right exists, put to stack
            pinin<-pinin+1
            pino[pinin]<-right[cur]
        }
        while (left[cur]>0){    #go to leaf and put right nodes to stack
             cur<-left[cur]
             if (left[cur]>0){  #if not leaf
                vecs<-allokoi.new(cur,vecs,lst,left,right,sibord)
             }
             if (right[cur]>0){ #if right exists, put to stack
                pinin<-pinin+1
                pino[pinin]<-right[cur]
             }
        }
    }
} 
      
if (dens) firstlevel<-0 else firstlevel<-min(lst$level)
xlim<-c(min(vecs[,1]),max(vecs[,2]))
ylim<-c(firstlevel,max(lst$level))
plot(x="",y="",xlab="",ylab="",xlim=xlim,ylim=ylim)

for (i in 1:itemnum){

    yc<-vecs[i,3]

    pare<-lst$parent[i]
    if (pare==0) lowlev<-firstlevel else lowlev<-lst$level[pare]

    segments(vecs[i,1],lowlev,vecs[i,1],yc)#,col=col,lwd=thick)
    segments(vecs[i,2],lowlev,vecs[i,2],yc)#,col=col,lwd=thick)

    if (left[i]==0){  #we are in leaf

       segments(vecs[i,1],yc,vecs[i,2],yc)#,col=col,lwd=thick)

    }
    else{

       childnum<-1
       curchi<-mt$child[i]
       while (mt$sibling[curchi]!=0){
           curchi<-mt$sibling[curchi]
           childnum<-childnum+1
       }

       sibpointer<-matrix(0,childnum,1)
       curchi<-mt$child[i]
       sibpointer[sibord[curchi]]<-curchi
       while (mt$sibling[curchi]!=0){
           curchi<-mt$sibling[curchi]
           sibpointer[sibord[curchi]]<-curchi
       }

       curchi<-sibpointer[1]
       x1<-vecs[curchi,1]      
       segments(vecs[i,1],yc,x1,yc)#,col=col,lwd=thick)
       x0<-vecs[curchi,2]

       cn<-2
       while (cn<=childnum){
             curchi<-sibpointer[cn]
             x1<-vecs[curchi,1] 
             segments(x0,yc,x1,yc)#,col=col,lwd=thick)
             x0<-vecs[curchi,2] 
             cn<-cn+1
       }

       segments(x0,yc,vecs[i,2],yc)#,col=col,lwd=thick)

    }
}


}





