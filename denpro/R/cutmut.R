cutmut<-function(mut,level,levels){
#
roots<-mut$roots
child<-mut$child
sibling<-mut$sibling
#
itemnum<-length(child)
rootnum<-length(roots)       
#
newroots<-matrix(0,itemnum,1)
newsibling<-sibling
ind<-0
#
for (i in 1:rootnum){
      curroot<-roots[i]
      pino<-matrix(0,itemnum,1)
      pino[1]<-curroot
      pinin<-1
      while (pinin>0){
          cur<-pino[pinin]      #take from stack
          pinin<-pinin-1
          # 
          # if cur acrosses the level, make cur root
          if (levels[cur]>level){
             ind<-ind+1
             newroots[ind]<-cur     # add to list
             newsibling[cur]<-0     # remove siblings
          }               
          # put to the stack
          if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
          }
          # go to left and put right nodes to the stack
          # if we have not already crossed the level
          while ((child[cur]>0) && (levels[cur]<=level)){
             cur<-child[cur]
             # if cur acrosses the level, make cur root
             if (levels[cur]>level){
                ind<-ind+1
                newroots[ind]<-cur  # add to list
                newsibling[cur]<-0     # remove siblings
             }  
             if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
             }
           }
       }
}                                        
if (ind>0){
    newroots<-newroots[1:ind]
}
else{
    newroots<-NULL
}
return(list(roots=newroots,sibling=newsibling))
}





