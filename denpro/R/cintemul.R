cintemul<-function(roots,child,sibling,volume,level){
#
#integrate function, over the level of roots, in the region of roots
#
itemnum<-length(child)
rootnum<-length(roots)
inte<-0
for (i in 1:rootnum){
    pino<-matrix(0,itemnum,1)
    valpino<-matrix(0,itemnum,1)  #level of parent
    pino[1]<-roots[i]
    valpino[1]<-0
    sibling[roots[i]]<-0
    #    
    pinin<-1
    while (pinin>0){
        cur<-pino[pinin]      #take from stack
        valcur<-valpino[pinin] 
        pinin<-pinin-1
        #
        if (level[cur]>0){
           inte<-inte+(level[cur]-max(valcur,0))*volume[cur]
        }
        #
        if (sibling[cur]>0){
              pinin<-pinin+1
              pino[pinin]<-sibling[cur]
              valpino[pinin]<-valcur
        }
        while (child[cur]>0){    #go to left and put right nodes to stack
              valcur<-level[cur]
              cur<-child[cur]
              #
              if (level[cur]>0){
                 inte<-inte+(level[cur]-max(valcur,0))*volume[cur]
              }
              #
              if (sibling[cur]>0){  #if cur has siblings
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
                 valpino[pinin]<-valcur
             }
        }
    }
}
#
return(inte)
}

