lowupp<-function(tr)
{
   nodenum<-length(tr$left)
   d<-length(tr$N)
   low<-matrix(0,nodenum,d)
   upp<-matrix(0,nodenum,d)

   upp[1,]<-tr$N  

   pino<-matrix(0,nodenum,1)
   pinoin<-1
   pino[pinoin]<-1
   while (pinoin>0){      # go through the nodes of tr2
       node<-pino[pinoin]
       pinoin<-pinoin-1

       while (tr$left[node]>0){   
       
          direk<-tr$direc[node]
          split<-tr$split[node]

          leftnode<-tr$left[node]
          rightnode<-tr$right[node]
  
          low[leftnode,]<-low[node,]
          upp[leftnode,]<-upp[node,]
          upp[leftnode,direk]<-tr$split[node]

          low[rightnode,]<-low[node,]
          low[rightnode,direk]<-tr$split[node]
          upp[rightnode,]<-upp[node,]

          # put right node to the stack (if exists)

          pinoin<-pinoin+1
          pino[pinoin]<-rightnode

          # go to left

          node<-leftnode
       }
   }

return(list(low=low,upp=upp))
}

