siborder.new<-function(mt)
{
#mt is multitree

itemnum<-length(mt$child)
sibord<-matrix(0,itemnum,1)

#order first roots

rootnum<-length(mt$roots)
for (i in 1:rootnum) sibord[mt$roots[i]]<-i

# then order the other

for (i in 1:rootnum){
   curroot<-mt$roots[i]
   if (mt$child[curroot]>0){
      pino<-matrix(0,itemnum,1)
      pino[1]<-mt$child[curroot]
      pinin<-1
      while (pinin>0){
          cur<-pino[pinin]      #take from stack
          pinin<-pinin-1
          # if not yet ordered, order siblings
          if (sibord[cur]==0){
              indu<-1
              sibord[cur]<-indu
              runner<-cur
              while (mt$sibling[runner]>0){
                  sibord[mt$sibling[runner]]<-indu
                  indu<-indu+1  
                  runner<-mt$sibling[runner] 
              }
          }
          # put to the stack 
          if (mt$sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-mt$sibling[cur]
          }
          # go to left and put right nodes to the stack
          while (mt$child[cur]>0){   
             cur<-mt$child[cur]
             # if not yet ordered, order siblings
             if (sibord[cur]==0){
                 indu<-1
                 sibord[cur]<-indu
                 runner<-cur
                 while (mt$sibling[runner]>0){
                     sibord[mt$sibling[runner]]<-indu
                     indu<-indu+1  
                     runner<-mt$sibling[runner]
                 }
             }
             if (mt$sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-mt$sibling[cur]
             }
           }
       }
   }
}                                  
return(sibord)
}















