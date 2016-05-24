final.table = function(wilcox=NULL, deseq=NULL, mageck=NULL, dataset, namecolumn=1, norm.function=median, type="genes", extractpattern = expression("^(.+?)_.+"))
{
  
  if(type=="genes")
  {
    # load dataset and walk through
    if(is.null(dataset$Gene))
    {data.frame$Gene  = as.character(sub(extractpattern,"\\1",data.frame[,namecolumn],perl=TRUE))}
    
    dataset.table = data.frame(
      name = dataset[,"Gene"],
      stringsAsFactors=FALSE)
    
    deseq = as.data.frame(deseq$genes)
    mageck = as.data.frame(mageck$genes)
    
    # load Wilcox info
    dataset.table$wilcox.FC = apply(dataset.table,1, function(y)
    {
      return(log2(wilcox[rownames(wilcox) == y["name"],"foldchange"]) )
    })
    
    dataset.table$wilcox.pval = apply(dataset.table,1, function(y)
    {
      return(wilcox[rownames(wilcox) == y["name"],"p.value"] )
    })
    
    # Load DESeq2 Info
    dataset.table$deseq2.FC = apply(dataset.table,1, function(y)
    {
      return(deseq[rownames(deseq) == y["name"],"log2FfoldChange"] )
    })
    
    dataset.table$deseq2.pval = apply(dataset.table,1, function(y)
    {
      return(deseq[rownames(deseq) == y["name"],"padj"] )
    })
    
    # Load MAGeCK info
      dataset.table$mageck.rank.pos = apply(dataset.table,1, function(y)
      {
        return(mageck[mageck$genes == y["name"],"rank.pos"] )
      })
      
      dataset.table$mageck.fdr.pos = apply(dataset.table,1, function(y)
      {
        return(mageck[mageck$genes == y["name"],"pos"] )
      })
      
      dataset.table$mageck.rank.neg = apply(dataset.table,1, function(y)
      {
        return(mageck[mageck$genes == y["name"],"rank.neg"] )
      })
      
      dataset.table$mageck.fdr.neg = apply(dataset.table,1, function(y)
      {
        return(mageck[mageck$genes == y["name"],"neg"] )
      })
    
    
  }
#   else if (type=="sgrna" || type=="all")
#   {
#     # load dataset and walk through
#     if(is.null(dataset$Gene))
#     {data.frame$Gene  = as.character(sub(extractpattern,"\\1",data.frame[,namecolumn],perl=TRUE))}
#     
#     dataset.table = data.frame(
#       name = dataset[,namecolumn],
#       stringsAsFactors=FALSE)
#     
#     deseq = as.data.frame(deseq$sgRNA)
#     mageck = as.data.frame(mageck$sgRNA)
#     
#     # use gene name information to load gene-based data
#     
#     # load Wilcox info
#     dataset.table$wilcox.FC = apply(dataset.table,1, function(y)
#     {
#       return(log2(wilcox[rownames(wilcox) == y["name"],"foldchange"]) )
#     })
#     
#     dataset.table$wilcox.pval = apply(dataset.table,1, function(y)
#     {
#       return(wilcox[rownames(wilcox) == y["name"],"p.value"] )
#     })
#     
#     # Load DESeq2 Info
#     dataset.table$deseq2.FC = apply(dataset.table,1, function(y)
#     {
#       return(deseq[deseq$genes == y["name"],"log2FfoldChange"] )
#     })
#     str(deseq)
#     dataset.table$deseq2.pval = apply(dataset.table,1, function(y)
#     {
#       return(deseq[deseq$genes == y["name"],"padj"] )
#     })
#     
#     # Load MAGeCK info
#     dataset.table$mageck.rank.pos = apply(dataset.table,1, function(y)
#     {
#       return(mageck[mageck$genes == y["name"],"rank.pos"] )
#     })
#     
#     dataset.table$mageck.fdr.pos = apply(dataset.table,1, function(y)
#     {
#       return(mageck[mageck$genes == y["name"],"pos"] )
#     })
#     
#     dataset.table$mageck.rank.neg = apply(dataset.table,1, function(y)
#     {
#       return(mageck[mageck$genes == y["name"],"rank.neg"] )
#     })
#     
#     dataset.table$mageck.fdr.neg = apply(dataset.table,1, function(y)
#     {
#       return(mageck[mageck$genes == y["name"],"neg"] )
#     })
#     
#   }
#   
  
  # final output genes
  # final output sgRNAs
  
  dataset.return=dataset.table
  

  # Final overall including genes and sgRNAs
  
  
  return(dataset.return)
  
}