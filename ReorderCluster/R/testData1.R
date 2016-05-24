testData1<-function()
{
  num_gene=400
  num_exp=100
  num=10
  var=30
  p=0.01
  mat=runif(num_gene*num_exp)
  dim(mat)=c(num_gene,num_exp)
  class=numeric(num_gene)
  label=c(1,2,3,4)
  for(i in 1:8)
  {
    start=1+50*(i-1)
    end=start+50-1
    s=1+(i-1)*num
    e=s+var-1
    
    t=runif((end-start+1)*(e-s+1))
    index=which(t<=p)
    t=rep(-1,(end-start+1)*(e-s+1))
    pp=runif(length(index))
    for(j in 1:length(index))
    {
      if(pp[j]<0.5) 
        t[index[j]]=0
      else
        t[index[j]]=1
    }
    mat[start:end,s:e]=t
    #if(i<=3)
    #class[start:end]=label[i]
    #else
    #{
    #  if(i<=6)
    #class[start:end]=label[i-3] 
    #  else
    #    class[start:end]=label[i-6] 
    #}
  }
  
  class[1:125]=1
  class[126:275]=2
  class[276:400]=3
 
  
  image(mat)
  browser()
  mat=cbind(mat,class)
  colnames(mat)=c(paste("F",1:num_exp,sep=""),"class")
  mat
}