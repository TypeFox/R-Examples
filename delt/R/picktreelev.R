picktreelev<-function(treeseq,leafnum){

tree<-treeseq$tree
#delnodes<-treeseq$delnodes
#delend<-treeseq$delend
leafs<-treeseq$leafs
indeksi<-detsi(leafs,leafnum)
endi<-treeseq$delnodeend[indeksi]
if (endi>0){       #if there is something to remove
  indeksit<-treeseq$delnodes[1:endi]
  re<-remnodes(tree$left,tree$right,indeksit)
  tree$left<-re$left
  tree$right<-re$right
}                                  
return(tree)
}
