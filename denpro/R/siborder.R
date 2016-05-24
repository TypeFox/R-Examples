siborder<-function(mt,crit,centers)
{
#mt is multitree

roots<-mt$roots
child<-mt$child
sibling<-mt$sibling

itemnum<-length(child)
sibord<-matrix(0,itemnum,1)

#order first roots

rootnum<-length(roots)
if (rootnum==1){
  sibord[roots[1]]<-1
}
else{
  rootlink<-matrix(0,itemnum,1)
  for (i in 1:(rootnum-1)){
     inde<-roots[i]
     rootlink[inde]<-roots[i+1]
  }
  sibord<-levord(roots[1],rootlink,sibord,centers,crit)
}

# then order the other

for (i in 1:rootnum){
   curroot<-roots[i]
   if (child[curroot]>0){
      pino<-matrix(0,itemnum,1)
      pino[1]<-child[curroot]
      pinin<-1
      while (pinin>0){
          cur<-pino[pinin]      #take from stack
          pinin<-pinin-1
          # if not yet ordered, order siblings
          if (sibord[cur]==0){
              sibord<-levord(cur,sibling,sibord,centers,crit)
          }
          # put to the stack 
          if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
          }
          # go to left and put right nodes to the stack
          while (child[cur]>0){   
             cur<-child[cur]
             # if not yet ordered, order siblings
             if (sibord[cur]==0){
                sibord<-levord(cur,sibling,sibord,centers,crit)
             }
             if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
             }
           }
       }
   }
}                                  
return(sibord)
}















