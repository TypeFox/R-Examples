levord<-function(beg,sibling,sibord,centers,crit){
#order at the given level
#
# find first
#
itemnum<-length(sibling)
diffe<-matrix(NA,itemnum,1)    #NA is infty
cur<-beg
curre<-centers[,cur]
diffe[cur]<-etais(curre,crit)
sibnum<-1   #if beg has no siblings, then sibnum=1 (beg is its own sibling)
while (sibling[cur]>0){
     cur<-sibling[cur] 
     curre<-centers[,cur]
     diffe[cur]<-etais(curre,crit)
     sibnum<-sibnum+1
}
first<-omaind(-diffe)
#sibord[first]<-1
#
# find distances to first
#
firstcenter<-centers[,first]
distofir<-matrix(NA,itemnum,1)
cur<-beg
curre<-centers[,cur]
distofir[cur]<-etais(curre,firstcenter)
while (sibling[cur]>0){
    cur<-sibling[cur]
    curre<-centers[,cur]
    distofir[cur]<-etais(curre,firstcenter)
}
#  
# fill sibord in the order of closest to first 
#
ind<-1
remain<-sibnum
while (remain>0){  
     nex<-omaind(distofir)
     sibord[nex]<-ind
     distofir[nex]<-NA 
     ind<-ind+1
     remain<-remain-1     
}
return(sibord)
}











