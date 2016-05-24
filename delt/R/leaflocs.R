leaflocs<-function(left,right){
#Finds location of leafs in binary tree, living in vector

itemnum<-length(left)
leafloc<-matrix(0,itemnum,1)
pino<-matrix(0,itemnum,1)
pino[1]<-1     #pino[1]<-root
pinin<-1
leafnum<-0
while (pinin>0){
    cur<-pino[pinin]      #take from stack
    pinin<-pinin-1
    if (left[cur]==0){    #if we are in leaf
       leafnum<-leafnum+1
       leafloc[leafnum]<-cur
    }
    else{
       while (left[cur]>0){  #go to leaf and put right nodes to stack
           pinin<-pinin+1
           pino[pinin]<-right[cur]
           cur<-left[cur]
       }
       #now we know we are in leaf
       leafnum<-leafnum+1
       leafloc[leafnum]<-cur  
    }
}
leafloc<-leafloc[1:leafnum]
return(list(leafloc=leafloc,leafnum=leafnum))
} 




