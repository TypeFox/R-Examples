touchstep.boundary<-function(node,curroot,boundrec,child,sibling,infopointer,
low,upp,rho=0)
{
# Checks whether "node" touches some of the leafs of the branch whose
# root is "curroot". Goes through the branch starting at "curroot".
# "comprec" is associated with the "node"
# "currec" is the bounding box of "cur"
# "pointrec" is associated with "cur"

d<-length(low[1,])
comprec<-matrix(0,2*d,1)
note<-infopointer[node]   #nodefinder[infopointer[node]]
for (i in 1:d){
    comprec[2*i-1]<-low[note,i]   #index[infopointer[node],i]
    comprec[2*i]<-upp[note,i]     #index[infopointer[node],i]
}

itemnum<-length(child)
pino<-matrix(0,itemnum,1)

potetouch<-1
istouch<-0
pino[1]<-curroot
pinin<-1
while ((pinin>0) && (istouch==0)){
      cur<-pino[pinin]      #take from stack
      pinin<-pinin-1

      # create currec and pointrec
      currec<-boundrec[cur,]
      pointrec<-matrix(0,2*d,1)
      note<-infopointer[cur]   #nodefinder[infopointer[cur]]
      for (i in 1:d){
         pointrec[2*i-1]<-low[note,i] #index[infopointer[cur],i]
         pointrec[2*i]<-upp[note,i]   #index[infopointer[cur],i]
      }
      # find touches                        
      potetouch<-touchi.boundary(comprec,currec,rho) 
      istouch<-touchi.boundary(comprec,pointrec,rho)

      # put to the stack
      if (sibling[cur]>0){
            pinin<-pinin+1
            pino[pinin]<-sibling[cur]
      }
      # go to left and put right nodes to the stack
      while ((child[cur]>0) && (istouch==0) && (potetouch==1)){
            cur<-child[cur]

            # create currec and pointrec
            currec<-boundrec[cur,]
            pointrec<-matrix(0,2*d,1)
            note<-infopointer[cur]   #nodefinder[infopointer[cur]]
            for (i in 1:d){
               pointrec[2*i-1]<-low[note,i] #index[infopointer[cur],i]
               pointrec[2*i]<-upp[note,i]   #index[infopointer[cur],i]
            }
            # find touches                        
            potetouch<-touchi.boundary(comprec,currec,rho) 
            istouch<-touchi.boundary(comprec,pointrec,rho)
 
            if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
            }
      }
}

return(istouch)
}                                    
