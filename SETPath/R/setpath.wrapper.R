setpath.wrapper <-
function(d1,d2,pathwaygenes,pathwaynames,M=1,transform=NULL,minalpha=NULL,normalize=TRUE,pvalue="chisq",npermutations=10000)
{
  K = length(pathwaynames)
  results = matrix(NA,K,2*(M+1)+2)
  dimnames(results)[[1]] = pathwaynames
  dimnames(results)[[2]] = c("n.genes",paste("alpha.0",1:M,sep="."),"T.0",paste("alpha.0",1:M,sep="."),"T.0","pval")
  
  # make sure gene names align:
  if(!identical(dimnames(d1)[[2]],dimnames(d2)[[2]]))
  {
    stop("d1 and d2 have different feature (column) names.")
  }
  # now run the method on each pathway:
  for(k in 1:K)
  {  
    # check formatting of pathways:
    missinggenes = setdiff(pathwaygenes[[k]],dimnames(d1)[[2]])
    if(length(missinggenes)>0)
    {
      warning(c("The following pathway genes are missing from the dataset:",missinggenes))
      pathwaygenes[[k]] = intersect(pathwaygenes[[k]],dimnames(d1)[[2]])
    }
    temp = setpath(d1[,pathwaygenes[[k]]],d2[,pathwaygenes[[k]]],M=M,transform=transform,verbose=TRUE,minalpha=minalpha,normalize=normalize,pvalue=pvalue,npermutations=npermutations)
    results[k,] = c(length(pathwaygenes[[k]]),temp$stats[,1],temp$stats[,2],temp$pval)
  }
  return(results)
}
