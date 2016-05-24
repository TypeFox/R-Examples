unmapped.genes = function(data, namecolumn=1, fullmatchcolumn=2, genes=NULL, extractpattern=expression("^(.+?)_.+")) 
{
  # input
  #   data = dataset sgRNA with readcount
  #   namecolumn
  #   fullmatchcolumn
  #   genes = NULL -> will list genes with unmapped sgRNAs
  #   genes = vector with gene names -> will list unmapped sgRNAs by name for the provided gene names
  
  
  # desired output
  # name of gene, unmapped sgRNAs
  # or
  # if gene name is provided
  # name of sgRNAs that are unmapped for this gene
  if(is.null(genes))
  {
    # get gene names of unmapped sgRNAs
    names = sub(extractpattern,"\\1",data[data[,fullmatchcolumn] == 0, namecolumn],perl=TRUE)
    if(length(names) > 0)
    {
      df.unmapped.sgRNA = data.frame(
        name = names,
        y = 1,
        stringsAsFactors=FALSE)
      
      # by aggregating, sum up how many sgRNAs for each gene were unmapped (readcount == 0)
      df.unmapped.sgRNA = aggregate(df.unmapped.sgRNA$y, by=list(df.unmapped.sgRNA$name), function(x) return(length(x[x==1])))
      colnames(df.unmapped.sgRNA) = c("name","sgRNA")
      rownames(df.unmapped.sgRNA) = df.unmapped.sgRNA$name
      
      return(df.unmapped.sgRNA)
    }
    else
    {
      df.unmapped.sgRNA = data.frame( name = "none", "sgRNA" = 0)
      colnames(df.unmapped.sgRNA) = c("name","sgRNA")
  
      
      return(df.unmapped.sgRNA)
    }
    
    
    
  }
  else
  {
    # get single sgRNA names for a given gene
    df.unmapped.sgRNA = data.frame(
      sgRNA = data[ data[,fullmatchcolumn] == 0 & sub(extractpattern,"\\1",data[, namecolumn],perl=TRUE) %in% genes, namecolumn],
      stringsAsFactors=FALSE)
    df.unmapped.sgRNA$gene = sub(extractpattern,"\\1",df.unmapped.sgRNA[, "sgRNA"],perl=TRUE) 
    
    return(df.unmapped.sgRNA)
  }

}