findleafs<-function(left,right)
{
# Finds location of leafs in binary tree, living in vector
# left, right are itemnum-vectors

# returns itemnum-vector, 1 in the location of nodes 0 non-terminal
# NA in positions not belonging to tree

# vector where binary tree is living may be larger than cardinality
# of nodes of the tree

itemnum<-length(left)
leafloc<-matrix(NA,itemnum,1)
pino<-matrix(0,itemnum,1)
pino[1]<-1     #pino[1]<-root
pinin<-1
while (pinin>0){
    cur<-pino[pinin]      #take from stack
    pinin<-pinin-1
    if (left[cur]==0){    #if we are in leaf
       leafloc[cur]<-1
    }
    else{
       while (left[cur]>0){  #go to leaf and put right nodes to stack
           leafloc[cur]<-0
           pinin<-pinin+1
           pino[pinin]<-right[cur]
           cur<-left[cur]
       }
       leafloc[cur]<-1  #now we know we are in leaf
    }
}
return(leafloc)
} 







