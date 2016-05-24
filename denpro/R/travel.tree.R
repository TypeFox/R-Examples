travel.tree<-function(parent,node)
{
mt<-multitree(parent) #roots<-mt$roots child<-mt$child sibling<-mt$sibling
itemnum<-length(parent)
nodes<-matrix(0,itemnum,1)

curroot<-node
counter<-0
if (mt$child[curroot]>0){
   pino<-matrix(0,itemnum,1)
   pino[1]<-mt$child[curroot]
   pinin<-1
   while (pinin>0){
        cur<-pino[pinin]      #take from stack
        pinin<-pinin-1
        counter<-counter+1
        nodes[counter]<-cur   
        # put to the stack 
        if (mt$sibling[cur]>0){
            pinin<-pinin+1
            pino[pinin]<-mt$sibling[cur]
        }
        # go to left and put right nodes to the stack
        while (mt$child[cur]>0){   
            cur<-mt$child[cur]
            counter<-counter+1
            nodes[counter]<-cur   
            if (mt$sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-mt$sibling[cur]
            }
        }
   }#while (pinin>0)
}
       
nodes<-nodes[1:counter]
return(nodes)
}

