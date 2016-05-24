colobary.nodes<-function(parent,nodes,paletti=NULL)
{
mt<-multitree(parent) #roots<-mt$roots child<-mt$child sibling<-mt$sibling

itemnum<-length(mt$child)
rootnum<-length(mt$roots)
if (is.null(paletti))
 paletti<-c("red","blue","green",
 "orange","navy","darkgreen",
 "orchid","aquamarine","turquoise",
 "pink","violet","magenta","chocolate","cyan",
 colors()[50:657],colors()[50:657])

colot<-matrix("",itemnum,1)
col<-colobary(parent,paletti)

nodenum<-length(parent)
allnodes<-matrix(0,nodenum,1)
#allnodes[1:length(nodes)]<-nodes
counter<-0  #length(nodes)
for (i in 1:length(nodes)){
    node<-nodes[i]
    tt<-travel.tree(parent,node)
    allnodes[(counter+1):(counter+length(tt))]<-tt
    counter<-counter+length(tt)
}
allnodes<-allnodes[1:counter]

for (i in 1:rootnum){
   curroot<-mt$roots[i]
   colot[curroot]<-col[curroot]  #"grey"  #paletti[hep]
   if (mt$child[curroot]>0){
      pino<-matrix(0,itemnum,1)
      pino[1]<-mt$child[curroot]
      pinin<-1
      while (pinin>0){
          cur<-pino[pinin]      #take from stack
          pinin<-pinin-1
          if (sum(cur==allnodes)>0)
              colot[cur]<-colot[parent[cur]]
          else{ 
                colot[cur]<-col[cur]  #"grey"  #paletti[hep]
          } 
          # put to the stack 
          if (mt$sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-mt$sibling[cur]
          }
          # go to left and put right nodes to the stack
          while (mt$child[cur]>0){   
             cur<-mt$child[cur]
             if (sum(cur==allnodes)>0)
                       colot[cur]<-colot[parent[cur]]
             else{ 
                    colot[cur]<-col[cur]  #"grey"  #paletti[hep]
             } 
             if (mt$sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-mt$sibling[cur]
             }
           }
       }#while (pinin>0)
   }
}       

ind<-setdiff(seq(1:itemnum),allnodes)
colot[ind]<-"grey"
                    
return(colot)
}


