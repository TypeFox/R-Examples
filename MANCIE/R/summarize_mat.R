summarize_mat<-function(mat_main,ann_main,mat_supp,ann_supp,n_limit=50,extend=100000,method="pca")
{
  mat_main=as.matrix(mat_main)
  mod_main=mat_main
  mat_supp=as.matrix(mat_supp)
  
  for (i in 1:dim(mat_main)[1]) 
  {
    ids=get_nearby(i,ann_supp,ann_main,n_limit=n_limit,extend=extend)
    if (length(ids)==0) {next}
    if (length(ids)==1)
    {
      PCA=mat_supp[ids,]
    }else
    {
      if (method=="max")
      {
        ids_max=ids[which.max(abs(cor(mat_main[i,],t(mat_supp[ids,]))))]
        PCA=mat_supp[ids_max,]
      }else
      {
        PCA=prcomp(t(mat_supp[ids,]),.scale=T)$x[,1]
      }
    }
    if (cor(PCA,mat_main[i,])<0) {PCA=-PCA}
    mod_main[i,]=PCA
  }
  
  mod_main
}
