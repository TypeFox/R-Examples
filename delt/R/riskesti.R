riskesti<-function(treeseq,n){
#Estimates risk for every alpha
#
#treeseq is list(tree,leafs,alfa,...)
#  tree is list(volum,nelem,...)
#n is the sample size
#
#Returns an alphalkm-vector
#
tr<-treeseq$tree
left<-tr$left
right<-tr$right
alfalkm<-length(treeseq$alfa)
toremove<-treeseq$delnodes
if (dim(t(toremove))[1]==1) maxrem<-1 else maxrem<-length(toremove[1,])
#mita jos toremove on skalaari?    
#
inum<-length(tr$vec)
ykk<-rep(1,inum)
nelem<-tr$nelem
volum<-tr$volum
#  estimated risk is sum of the info over leafs, info is vector which
#  we have to sum over leafs 
info<-nelem*(ykk-nelem*(1+1/n))/volum  #/n^2
#
risks<-matrix(0,alfalkm,1)
risks[1]<-leafsum(info,root=1,left,right) #kun alpha=0, ei ole poist mitaan
cursum<-risks[1]
for (i in 1:alfalkm){
    if (maxrem==1){ 
       poista<-toremove[i]
       sumsubtree<-leafsum(info,root=poista,left,right)
       cursum<-cursum-sumsubtree+info[poista]
       left[poista]<-0
       right[poista]<-0
    }
    else{
      j<-1
      while ((j<=maxrem) && (toremove[i,j]>0)){
         poista<-toremove[i,j]
         sumsubtree<-leafsum(info,root=poista,left,right)
         cursum<-cursum-sumsubtree+info[poista]
         left[poista]<-0
         right[poista]<-0  
         j<-j+1
      }
    }
    risks[i]<-cursum
}
return(t(risks))
}




