eval.pick<-function(treeseq,leaf){
#
tree<-treeseq$tree
delnodes<-treeseq$delnodes
delend<-treeseq$delend
leafs<-treeseq$leafs
indeksi<-detsi(leafs,leaf)
endi<-delend[indeksi]
if (endi>0){       #if there is something to remove
  indeksit<-delnodes[1:endi]
  re<-remnodes(tree$left,tree$right,indeksit)
  tree$left<-re$left
  tree$right<-re$right
}  

####################################################
cl<-length(tree$left)
d<-length(tree$N)

down<-matrix(0,cl,d)
high<-matrix(0,cl,d)
for (i in 1:cl){
   for (j in 1:d){
      down[i,j]<-tree$low[i,j]+1
      high[i,j]<-tree$upp[i,j]
    }
}

ll<-leaflocs(tree$left[1:cl],tree$right[1:cl])
leafloc<-ll$leafloc
leafnum<-ll$leafnum

value<-matrix(0,leafnum,1)

efek<-0
i<-1
while (i<=leafnum){  
   node<-leafloc[i]

   if (tree$mean[node]>0){
     efek<-efek+1

     value[efek]<-tree$mean[node]
 
     for (j in 1:d){
         down[efek,j]<-tree$low[node,j]
         high[efek,j]<-tree$upp[node,j]
     }
   }
   i<-i+1
}
value<-value[1:efek]
if (efek>1){
   down<-down[1:efek,]
   high<-high[1:efek,]
}
else{
   apudown<-matrix(0,1,d)
   apuhigh<-matrix(0,1,d)
   for (ddd in 1:d){
        apudown[1,ddd]<-down[1,ddd]
        apuhigh[1,ddd]<-high[1,ddd]
   }
   down<-apudown
   high<-apuhigh
}

###################################################
tree$value<-value
tree$down<-down
tree$high<-high
   
return(tree)
}
