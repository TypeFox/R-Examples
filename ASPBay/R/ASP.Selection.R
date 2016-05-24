ASP.Selection <-
function(Tem_Gen, Index_Gen, IBD,k=log(10000))
{
  if ( ncol(Tem_Gen)!=ncol(Index_Gen) ) return("Error : Controls and index cases don't have the same number of genotypes")
  if ( length(IBD)!=nrow(Index_Gen) ) return("Error : We need the same number of index cases and IBD states")
  
  p <- ncol(Tem_Gen)
  score <- ASP.Score(Tem_Gen, Index_Gen, IBD)$Value
  
  stat <- max(score**2) - score**2
  
  names(score) <- colnames(Tem_Gen)
  names(stat) <- colnames(Tem_Gen)
  
  return(list(score=score, dis_stat=stat, SNP_subset = (1:p)[stat<k], SNPnames_subset = colnames(Tem_Gen)[stat<k]))
}

