funMerge<-function(ind,row,col,hc,node,maxI,maxJ,cpp)
{
  
  if(!row%in%node[[ind]]$left)
  {
    hc$merge[ind,]=hc$merge[ind,c(2,1)]
  }
  if(cpp)
  {
    i1=maxI[row,col]+1
    j1=maxJ[row,col]+1
  }
  else
  {
    i1=maxI[row,col]
    j1=maxJ[row,col]
  }
  if(i1!=row)
  {
    hc=funMerge(hc$merge[ind,1],row,i1,hc,node,maxI,maxJ,cpp)
  }
  if(j1!=col)
  {
    hc=funMerge(hc$merge[ind,2],j1,col,hc,node,maxI,maxJ,cpp)
  }
  return(hc)	

}