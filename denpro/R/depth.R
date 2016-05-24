depth<-function(mt){
#finds the dephts of the nodes
#
#mt is a result from multitree or pruneprof
#
roots<-mt$roots
child<-mt$child
sibling<-mt$sibling
#
itemnum<-length(child)
depths<-matrix(0,itemnum,1)
#
rootnum<-length(roots)
#
for (i in 1:rootnum){
    pino<-matrix(0,itemnum,1)
    pino[1]<-roots[i]  
    pinin<-1
    depths[roots[i]]<-1
    while (pinin>0){
        cur<-pino[pinin]      #take from stack
        pinin<-pinin-1
        if (sibling[cur]>0){    #put right to stack
             pinin<-pinin+1
             pino[pinin]<-sibling[cur]
             depths[sibling[cur]]<-depths[cur]
        }
        while (child[cur]>0){    #go to leaf and put right nodes to stack
             chi<-child[cur]
             depths[chi]<-depths[cur]+1
             if (sibling[chi]>0){ 
                   pinin<-pinin+1
                   pino[pinin]<-sibling[chi]
                   depths[sibling[chi]]<-depths[cur]+1
             }
             cur<-chi
        }
    }
}
return(depths)
}







