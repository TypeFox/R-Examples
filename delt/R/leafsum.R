leafsum<-function(info,root,left,right){
#Calculates the sum of info over leaves of the subtree
#
#info is itemnum-vector
#root is in 1:itemnum
#left, right are itemnum-vectors, links to childs, 0 if leaf
#
itemnum<-length(info)
sum<-0
pino<-matrix(0,itemnum,1)
pino[1]<-root
pinin<-1
while (pinin>0){
  cur<-pino[pinin]
  pinin<-pinin-1
  if (left[cur]==0){  #if we are in leaf, sum
    sum<-sum+info[cur]
  }
  else{  
    while (left[cur]>0){  #put right on the stack and go to left 
       pinin<-pinin+1
       pino[pinin]<-right[cur]
       cur<-left[cur]
    }
    sum<-sum+info[cur]   #now we know we are in leaf      
  }
}
return(sum)
} 
