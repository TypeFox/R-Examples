pruneprof<-function(mt){
#prunes profile so that only root and nodes with siblings are left
#
#mt is a result from multitree
#
roots<-mt$roots
child<-mt$child
sibling<-mt$sibling
siborder<-mt$siborder
#
itemnum<-length(child)
newchild<-matrix(0,itemnum,1)
#
rootnum<-length(roots)
#
for (i in 1:rootnum){
    pino<-matrix(0,itemnum,1)
    pino[1]<-roots[i]  
    pinin<-1
    while (pinin>0){
        cur<-pino[pinin]      #take from stack
        pinin<-pinin-1
        if (sibling[cur]>0){
              pinin<-pinin+1
              pino[pinin]<-sibling[cur]
        }
        while (child[cur]>0){    #go to left and put right nodes to stack
             candi<-child[cur]
             while ((child[candi]>0) && (sibling[candi]==0)){
                 candi<-child[candi]
             }
             if (sibling[candi]>0){  #if candi has siblings
                newchild[cur]<-candi
                pinin<-pinin+1
                pino[pinin]<-sibling[candi]
             } 
             cur<-candi
        }
    }
}
return(list(roots=roots,child=newchild,sibling=sibling,siborder=siborder))
}







