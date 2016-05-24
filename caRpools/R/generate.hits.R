generate.hits = function(wilcox=NULL, deseq=NULL, mageck=NULL,  type="enriched", cutoff.deseq = 0.001, cutoff.wilcox = 0.05, cutoff.mageck = 0.05, cutoff.override=FALSE, cutoff.hits=NULL, plot.genes="overlapping")
{
  
  if(type=="enriched")
  {
    
    # separate by enrichment 
    if(!is.null(wilcox))
    {
      # wilcox
      df.wilcox = wilcox[wilcox$foldchange >= 1 ,c("foldchange","p.value")]
      
      df.wilcox = df.wilcox[order(
        df.wilcox$p.value,
        decreasing = FALSE 
      ),]
    }
    
    if(!is.null(deseq))
    {
      
      #DESeq2
      # Deseq make sure it is -1 * foldchange, as deseq compare UNTREATED vs TREATED< all other treated vs untreated!
      df.deseq = as.data.frame(deseq$genes[deseq$genes$log2FoldChange >= 0 ,c("log2FoldChange","padj")])
      
      df.deseq = df.deseq[order(
        df.deseq$padj,
        decreasing = FALSE
      ),]
    }
    
    
    if(!is.null(mageck))
    {
      # MAGeCK
      df.mageck = mageck$genes[, c("pos","rank.pos","sgrna.pos.good")]
      
      colnames(df.mageck) = c("fdr","rank","sgrna")
      df.mageck = df.mageck[order(
        df.mageck$rank,
        decreasing=FALSE,
        na.last=TRUE
      ),]
      
    }
    
    # Loading of data complete
    
    # Combine datasets using the cutoff
    # put in the rank of that gene in the data analysis table, if not present, put NA
    
    
    
  } else if(type=="depleted")
  {
    
    if(!is.null(wilcox))
    {
      
      # wilcox
      df.wilcox = wilcox[wilcox$foldchange < 1 , c("foldchange","p.value")]
      df.wilcox = df.wilcox[order(
        df.wilcox$p.value,
        decreasing = FALSE
      ),]
    }
    # print(df.wilcox)
    if(!is.null(deseq))
    {
      
      #DESeq2
      df.deseq = as.data.frame(deseq$genes[deseq$genes$log2FoldChange < 0 ,c("log2FoldChange","padj")])
      
      df.deseq = df.deseq[order(
        df.deseq$padj,
        decreasing = FALSE
      ),]
    }
    
    
    if(!is.null(mageck))
    {
      # MAGeCK
      df.mageck = mageck$genes[, c("neg","rank.neg","sgrna.neg.good")]
      
      colnames(df.mageck) = c("fdr","rank","sgrna")
      
      df.mageck = df.mageck[order(
        df.mageck$rank,
        decreasing=FALSE,
        na.last=TRUE
      ),]
      
    }
    
  } else {
    stop("No type selected. Please select enriched or depleted.")
  }
  #### Data loading complete
  
  # add ranking just by length of dataset
  #if(!is.null(riger)) {df.riger$rank = c(1:nrow(df.riger))}
  if(!is.null(wilcox)) {df.wilcox$rank = c(1:nrow(df.wilcox))}
  if(!is.null(deseq)) {df.deseq$rank = c(1:nrow(df.deseq))}
  
  # create joining data.frame that has all gene names, use any of the above dataset
  
  if(!is.null(wilcox)) { v.join = rownames(wilcox)
                           
  } else if(!is.null(deseq)) { v.join = rownames(deseq) 
                               
  } else if(!is.null(mageck)) { v.join = rownames(mageck) 
                                
  } else {stop("No datasets provided")}
  
  # Create output data.frame
  # rownames = gene names
  # then
  # EITHER p.value, foldchange of RIGER, wilcox, DESEQ + AUC of DSS
  # OR rank of RIGER, wilcox, DESEQ, DSS that are within cutoff
  
  # remove all genes that have NA on all of them
  # remove columns that only have NAs (means not data was provided)
  
  
  ## Prepare output of a list used for Readcount scatter with an output of gene names only of ALL analysis methods
  
  
  # introduce output variable
  df.output = NA
  
  # wilcox
  if(!is.null(wilcox)) {
    
    if(identical(cutoff.override,TRUE) && is.numeric(cutoff.hits))
    {
      df.wilcox = df.wilcox[1:cutoff.hits,]
    }
    else if(identical(cutoff.override,TRUE) && is.null(cutoff.hits))
    {
      df.wilcox = df.wilcox
    }
    else
    {
      df.wilcox = df.wilcox[df.wilcox$p.value <= cutoff.wilcox,]
    }
    
    #print(df.wilcox)
    
    #print(length(row(df.wilcox)))
    if(length(row(df.wilcox)) >=1)
    {
      for(i in 1:nrow(df.wilcox))
      {
        if(match(rownames(df.wilcox[i,]), df.output, nomatch = 1) == 1)
        {
          df.output = c(df.output,rownames(df.wilcox[i,]))
        }
        
      }
    }
    
  }
  
  
  
  # DESEQ2
  if(!is.null(deseq)) {
    
    if(identical(cutoff.override,TRUE) && is.numeric(cutoff.hits))
    {
      df.deseq = df.deseq[1:cutoff.hits,]
    }
    else if(identical(cutoff.override,TRUE) && is.null(cutoff.hits))
    {
      df.deseq = df.deseq
    }
    else
    {
      df.deseq = df.deseq[df.deseq$padj <= cutoff.deseq,]
    }
    
    #print(df.deseq)
    if(length(row(df.deseq) >= 1))
    {
      for(i in 1:nrow(df.deseq))
      {
        if(match(rownames(df.deseq[i,]), df.output, nomatch = 1) == 1)
        {
          df.output = c(df.output,rownames(df.deseq[i,]))
        }
        
      }
    }
    
  }
  
  # MAGeCK
  if(!is.null(mageck)) {
    
    if(identical(cutoff.override,TRUE) && is.numeric(cutoff.hits))
    {
      df.mageck = df.mageck[1:cutoff.hits,]
    }
    else if(identical(cutoff.override,TRUE) && is.null(cutoff.hits))
    {
      df.mageck = df.mageck
    }
    else
    {
      df.mageck = df.mageck[df.mageck$fdr <= cutoff.mageck,]
    }
    
    
    #print(df.mageck)
    if(length(row(df.mageck)) >= 1)
    {
      for(i in 1:nrow(df.mageck))
      {
        if(match(rownames(df.mageck[i,]), df.output, nomatch = 1) == 1)
        {
          df.output = c(df.output,rownames(df.mageck[i,]))
        }
        
      }
    }
  }
  
  # print(df.output)
  
  # by default: only overlapping hits will be plotted
  if(plot.genes=="overlapping") 
  {
    # Remove genes that are not present in all methods, so that only overlapping ones are plotted
    df.output = NULL
    # wilcox
    if(length(row(df.mageck)) >= 1 && length(row(df.wilcox)) >= 1 && length(row(df.deseq)) >= 1)
    {
      # check if hits are in other datasets
      for (i in 1:nrow(df.mageck))
      {
        
        #if(match(rownames(df.mageck[i,]), rownames(df.wilcox), nomatch = 0) == 1 && match(rownames(df.mageck[i,]), rownames(df.deseq), nomatch = 0) == 1)
        if(rownames(df.mageck[i,]) %in% rownames(df.deseq) && rownames(df.mageck[i,]) %in% rownames(df.wilcox))
        {
          
          df.output = c(df.output,rownames(df.mageck[i,]))
        }
      }
    }
    #else if(length(row(df.mageck)) >= 1) # Only plot Mageck
    #{
    #  df.output = rownames(df.mageck)
    #}
    
  }
  
  
  
  
  return(df.output)
  
  
}