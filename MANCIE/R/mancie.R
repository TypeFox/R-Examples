mancie<-function(mat_main,mat_supp,cutoff1=0.5,cutoff2=0)
{
  mat_main=as.matrix(mat_main)
  mat_supp=as.matrix(mat_supp)
  mat_final=mat_main
  cat1=cat2=cat3=0
  
  for (i in 1:dim(mat_main)[1]) # expression values are not changed if normalization fails
  {
    v1=mat_main[i,]
    v2=mat_supp[i,]
    
    if (cor(v1,v2)>cutoff1)
    {
      PCA=prcomp(cbind(v1,v2),.scale=T)$x[,1]
      if (cor(PCA,v1)<0) {PCA=-PCA}
      cat1=cat1+1
    }else if (cor(v1,v2)>cutoff2)
    {
      PCA=(scale(v1)+cor(v1,v2)*scale(v2))[,1]
      cat2=cat2+1
    }else
    {
      PCA=v1
      cat3=cat3+1
    }
    
    mat_final[i,]=(PCA-mean(PCA))/sd(PCA)*sd(v1)+mean(v1)
  }
  
  cat1=cat1/dim(mat_final)[1]
  cat2=cat2/dim(mat_final)[1]
  cat3=cat3/dim(mat_final)[1]
  cat(paste("Percentage of rows with correlation >cutoff1:",cat1,"\n"))
  if (cat1<0.1 || cat1>0.9) {warning("Too few or too many rows in this category!")}
  cat(paste("Percentage of rows with correlation <=cutoff1 and >cutoff2:",cat2,"\n"))
  if (cat2<0.1 || cat2>0.9) {warning("Too few or too many rows in this category!")}
  cat(paste("Percentage of rows with correlation <=cutoff2:",cat3,"\n"))
 
  mat_final
}
