dend2parent<-function(hc,dendat)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]
nodenum<-length(hc$height)+n
parent<-matrix(0,nodenum,1)
volume<-matrix(0,nodenum,1)
center<-matrix(0,d,nodenum)
level<-matrix(0,nodenum,1)
level[(n+1):nodenum]<-hc$height
pointers<-matrix(0,n,1)      #pointers dendat to tree 

joinnum<-length(hc$height)
lkm<-matrix(0,joinnum,1)
vec1<-matrix(0,1,d)
vec2<-matrix(0,1,d)
v1<-0
v2<-0
vapaa<-1
for (i in 1:joinnum){
   node1<-hc$merge[i,1]
   node2<-hc$merge[i,2]
   if (node1<0){
        parent[vapaa]<-i+n
        for (j in 1:d) vec1[j]<-dendat[-node1,j]
        v1<-1
        center[,vapaa]<-vec1
        volume[vapaa]<-v1
        pointers[-node1]<-vapaa
        lkm1<-1
        vapaa<-vapaa+1
   }
   else{
        parent[node1+n]<-i+n
        vec1<-center[,node1+n]
        v1<-volume[node1+n]
        lkm1<-lkm[node1]
   }
   if (node2<0){
        parent[vapaa]<-i+n
        for (j in 1:d) vec2[j]<-dendat[-node2,j]
        v2<-1
        center[,vapaa]<-vec2
        volume[vapaa]<-v2
        pointers[-node2]<-vapaa
        lkm2<-1
        vapaa<-vapaa+1
   }
   else{
        parent[node2+n]<-i+n
        vec2<-center[,node2+n]
        v2<-volume[node2+n]
        lkm2<-lkm[node2]
   }
   volume[i+n]<-1.1*(v1+v2)
   center[,i+n]<-(lkm1*vec1+lkm2*vec2)/(lkm1+lkm2)
   lkm[i]<-lkm1+lkm2
}

apoin<-matrix(0,n,1)
for (i in 1:n) apoin[i]<-nodenum-pointers[i]+1
apar<-parent[nodenum:1]
apar2<-matrix(0,n,1)
for (i in 1:nodenum) if (apar[i]!=0) apar2[i]<-nodenum-apar[i]+1

return(list(parent=apar2,level=level[nodenum:1],
volume=volume[nodenum:1],center=center[,nodenum:1],pointers=apoin))
}


