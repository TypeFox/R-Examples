allokoi.new<-function(cur,vecs,lst,left,right,sibord)
{
# allocates space for all children of "cur"

# Calculate the number of childs and sum of volumes of childs
now<-left[cur]
childnum<-1
childvolume<-lst$volume[now]
while (right[now]>0){
  now<-right[now]
  childnum<-childnum+1
  childvolume<-childvolume+lst$volume[now]
}
 
gaplen<-(lst$volume[cur]-childvolume)/(childnum+1)

if (childnum==1){
   now<-left[cur]
   xbeg<-gaplen+vecs[cur,1]
   xend<-xbeg+lst$volume[now]
   ycoo<-lst$level[now]
   vecs[now,]<-c(xbeg,xend,ycoo)
}
else{
  siblinks<-matrix(0,childnum,1)  #make siblinks in right order
  now<-left[cur] 
  sior<-sibord[now]
  siblinks[sior]<-now
  while (right[now]>0){
    now<-right[now]
    sior<-sibord[now]
    siblinks[sior]<-now
  }
  xend<-vecs[cur,1]      #initialize xend 
  for (i in 1:childnum){
     now<-siblinks[i]
     xbeg<-gaplen+xend
     xend<-xbeg+lst$volume[now]
     ycoo<-lst$level[now]
     vecs[now,]<-c(xbeg,xend,ycoo)
  }
} 


return(vecs)
}


