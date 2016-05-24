carpools.hit.scatter = function(wilcox=NULL, deseq=NULL, mageck=NULL, dataset, dataset.names = NULL, namecolumn=1, fullmatchcolumn=2, title="Read Count", xlab="Readcount Dataset1", ylab="Readcount Dataset2", labelgenes=NULL, labelcolor="orange", extractpattern=expression("^(.+?)_.+"), plotline=TRUE, normalize=TRUE, norm.function=median, offsetplot=1.2, center=FALSE, aggregated=FALSE, type="enriched", cutoff.deseq = 0.001, cutoff.wilcox = 0.05, cutoff.mageck = 0.05, cutoff.override=FALSE, cutoff.hits=NULL, plot.genes="overlapping", pch=16, col = rgb(0, 0, 0, alpha = 0.65))
{
  if(is.null(labelgenes))
  {
    
    # Call function to gene hit list
    
    df.output = generate.hits(wilcox=wilcox, deseq=deseq, mageck=mageck, type=type, cutoff.deseq = cutoff.deseq, cutoff.wilcox = cutoff.wilcox, cutoff.mageck = cutoff.mageck, cutoff.override=cutoff.override, cutoff.hits=cutoff.hits, plot.genes=plot.genes)

    # Call plot.read.count.vs to create paired plots
    #print(df.output)
    if(length(df.output) >= 1)
    {
      
      # Sorting 
      df.output = sort(df.output[!is.na(df.output)], na.last=TRUE, decreasing=FALSE) 
      
      for(i in 1:length(df.output))
      {
        carpools.read.count.vs(dataset=dataset, dataset.names = dataset.names, pairs=TRUE, namecolumn=namecolumn, fullmatchcolumn=fullmatchcolumn, title=title, pch=16, col = col, plotline=plotline, normalize=normalize, norm.function=norm.function, labelgenes=df.output[i], labelcolor=labelcolor, center=center, aggregated=aggregated, type=type)
      }
      
      return(df.output)
    }
    
    ###### End of creating automatic plots for all hits if labelgenes == NULL
  }
  
  else{
    
    ### Labelgene provided -> plot pairs for this gene
    # Call plot.read.count.vs to create paired plots
    carpools.read.count.vs(dataset=dataset, dataset.names = dataset.names, pairs=TRUE, namecolumn=namecolumn, fullmatchcolumn=fullmatchcolumn, title=title, pch=16, col = col, plotline=plotline, normalize=normalize, norm.function=norm.function, labelgenes=labelgenes, labelcolor=labelcolor, center=center, aggregated=aggregated, type=type)
    
  }
  
    
  
}