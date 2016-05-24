addnode<-function(inde,curre,curdep,left,right,parent,low,upp,N,numnode){
#(inde,curre,curdep,left,right,deplink,low,upp,enofatdep,N,numnode){
#
#inde is d-vector: index (gridpoint) to be added
#curre is pointer to vectors left,right,...
#
d<-length(inde)
apu<-depth2com(curdep,N)
curdir<-apu$direc
depatd<-apu$depind
depit<-log(N,base=2)
#depit[d]<-depit[d]+1
#
while (curdir<=(d-1)){
    ind<-inde[curdir]
    while (depatd<=depit[curdir]){
        mid<-(low[curre]+upp[curre])/2
        if (ind<=mid){
           left[curre]<-numnode+1
           parent[numnode+1]<-curre
           low[numnode+1]<-low[curre]
           upp[numnode+1]<-floor(mid)
        }
        else{
           right[curre]<-numnode+1
           parent[numnode+1]<-curre
           low[numnode+1]<-ceiling(mid)
           upp[numnode+1]<-upp[curre]
        }
        numnode<-numnode+1
        curre<-numnode
        depatd<-depatd+1
        curdep<-curdep+1
#        deplink[endofatdep[curdep]]<-numnode
#        deplink[numnode]<-0
#        endofatdep[curdep]<-numnode
    }
    #
    # Last node of this dimension (first node of next dimension)
    #
    curdir<-curdir+1
    ind<-inde[curdir]
    low[curre]<-1
    upp[curre]<-N[curdir]
    mid<-(low[curre]+upp[curre])/2
    if (ind<=mid){
           left[curre]<-numnode+1
           parent[numnode+1]<-curre
           low[numnode+1]<-low[curre]
           upp[numnode+1]<-floor(mid)
    }
    else{
           right[curre]<-numnode+1
           parent[numnode+1]<-curre
           low[numnode+1]<-ceiling(mid)
           upp[numnode+1]<-upp[curre]
    }
    depatd<-2
    numnode<-numnode+1
    curre<-numnode
    curdep<-curdep+1
#    deplink[endofatdep[curdep]]<-numnode
#    deplink[curre]<-0
#    endofatdep[curdep]<-numnode
}
#
# Last dimension 
#
ind<-inde[curdir]
while (depatd<=depit[curdir]){
        mid<-(low[curre]+upp[curre])/2
        if (ind<=mid){
           left[curre]<-numnode+1
           parent[numnode+1]<-curre
           low[numnode+1]<-low[curre]
           upp[numnode+1]<-floor(mid)
        }
        else{
           right[curre]<-numnode+1
           parent[numnode+1]<-curre
           low[numnode+1]<-ceiling(mid)
           upp[numnode+1]<-upp[curre]
        }
        numnode<-numnode+1
        curre<-numnode
        depatd<-depatd+1
        curdep<-curdep+1
#        deplink[endofatdep[curdep]]<-numnode
#        deplink[curre]<-0
#        endofatdep[curdep]<-numnode
}
#
# Last node of last dimension
#
#left[curre]<-0
#right[curre]<-0
#
#return(list(numnode=numnode,left=left,right=right,deplink=deplink,low=low,
#upp=upp,endofatdep=endofatdep))
return(list(numnode=numnode,left=left,right=right,parent=parent,low=low,
upp=upp,nodeloc=numnode))
}










