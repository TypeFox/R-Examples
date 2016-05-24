allokoi<-function(vecs,cur,child,sibling,sibord,levels,volumes)
{
#Finds coordinates of a node
#sibord,levels,volumes are nodenum-vector

# Calculate the number of childs and sum of volumes of childs
now<-child[cur]
childnum<-1
childvolume<-volumes[now]
while (sibling[now]>0){
  now<-sibling[now]
  childnum<-childnum+1
  childvolume<-childvolume+volumes[now]
}
 
gaplen<-(volumes[cur]-childvolume)/(childnum+1)

if (childnum==1){
   now<-child[cur]
   xbeg<-gaplen+vecs[cur,1]
   xend<-xbeg+volumes[now]
   ycoo<-levels[now]
   vecs[now,]<-c(xbeg,ycoo,xend,ycoo)
}
else{
  siblinks<-matrix(0,childnum,1)  #make siblinks in right order
  now<-child[cur] 
  sior<-sibord[now]
  siblinks[sior]<-now
  while (sibling[now]>0){
    now<-sibling[now]
    sior<-sibord[now]
    siblinks[sior]<-now
  }
  xend<-vecs[cur,1]      #initialize xend 
  for (i in 1:childnum){
     now<-siblinks[i]
     xbeg<-gaplen+xend
     xend<-xbeg+volumes[now]
     ycoo<-levels[now]
     vecs[now,]<-c(xbeg,ycoo,xend,ycoo)
  }
} 
return(vecs)
}



