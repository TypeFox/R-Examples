childcode<-function(level,runnumb){
#
#input: level and ordinal number in the listing of the full binary tree
#output: ordinal number of the left and right child
#
nodeNum<-0
i<-0
while (i<=(level-2)){
   nodeNum<-nodeNum+2^i
   i<-i+1
}
levBeg<-nodeNum+1
whichInLevel<-runnumb-levBeg+1
levBegNex<-levBeg+2^(level-1)
#
left<-levBegNex+(whichInLevel-1)*2
right<-left+1
#
return(list(left=left,right=right))
}
