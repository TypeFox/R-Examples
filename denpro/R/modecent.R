modecent<-function(lst){
#
parents<-lst$parent
levels<-lst$level
volumes<-lst$volume   
centers<-lst$center
d<-dim(centers)[1]            #d<-length(centers[,1])
#
mlkm<-moodilkm(parents)
modloc<-mlkm$modloc
lkm<-mlkm$lkm
#
mut<-multitree(parents)
roots<-mut$roots
child<-mut$child
sibling<-mut$sibling 
#
crit<-rep(0,d)               #order so that 1st closest to origo
sibord<-siborder(mut,crit,centers)   
#
itemnum<-length(parents)
vecs<-matrix(NA,itemnum,4)
vecs<-alloroot(vecs,roots,sibord,levels,volumes)
vecs<-plotdata(roots,child,sibling,sibord,levels,volumes,vecs) 
#
res<-matrix(0,lkm,d)
#
for (i in 1:lkm){
   sija<-modloc[i]
   res[i,]<-centers[,sija]
}
#
ord<-vecs[,1]   #in this order we want modes
ordpick<-matrix(0,lkm,1)
for (i in 1:lkm){
  sija<-modloc[i]
  ordpick[i]<-ord[sija]
}
#
pointer<-seq(1:lkm)
pointer<-omaord2(pointer,ordpick) #pointer on the order determined by ord
#
endres<-res
for (i in 1:lkm){
  sija<-pointer[i]
  endres[i,]<-res[sija,]
}
#
return(endres)
}

