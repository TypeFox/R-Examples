alloroot<-function(vecs,roots,sibord,levels,volumes)
{
rootnum<-length(roots)

# Calculate sum of volumes of roots
rootsvolume<-0
for (i in 1:rootnum){
  now<-roots[i]
  rootsvolume<-rootsvolume+volumes[now]
}

basis<-rootsvolume+rootsvolume/4
 
gaplen<-(basis-rootsvolume)/(rootnum+1)

rootlinks<-matrix(0,rootnum,1)  #make links in right order

if (rootnum==1) rootlinks[1]<-roots[1]  #1
else{
for (i in 1:rootnum){
  now<-roots[i]
  roor<-sibord[now]
  rootlinks[roor]<-now
}
}
xbeg<-0
xend<-0
for (i in 1:rootnum){
  now<-rootlinks[i]
  xbeg<-gaplen+xend
  xend<-xbeg+volumes[now]
  ycoo<-levels[now]
  vecs[now,]<-c(xbeg,ycoo,xend,ycoo)
}
return(vecs)
}
