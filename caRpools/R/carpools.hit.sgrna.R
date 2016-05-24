carpools.hit.sgrna = function(wilcox=NULL, deseq=NULL, mageck=NULL, dataset=NULL, dataset.names = NULL, namecolumn=1, fullmatchcolumn=2, norm.function=median, extractpattern=expression("^(.+?)_.+"), put.names=TRUE, type="enriched", labelgenes=NULL, cutoff.deseq = 0.05, cutoff.wilcox = 0.05, cutoff.mageck = 0.05, cutoff.override=FALSE, plot.genes="overlapping", cutoff.hits=NULL, plot.type=NULL, controls.target=NULL, controls.nontarget=NULL)
{
  
  
  # Call function to gene hit list
  if(is.null(labelgenes))
  {
    df.output = generate.hits(wilcox=wilcox, deseq=deseq, mageck=mageck, type=type, cutoff.deseq = cutoff.deseq, cutoff.wilcox = cutoff.wilcox, cutoff.mageck = cutoff.mageck, cutoff.override=cutoff.override, cutoff.hits=cutoff.hits, plot.genes=plot.genes)
    
  } else {
    df.output = labelgenes
  }
  
    #print(df.output)
    # Call plot.read.count.vs to create paired plots
    
    par(mfrow=c(1,1))
    par.old <- par
  
    if(is.null(plot.type))
    {
      par(mfrow=c(1,1))
    }
  
  
  if(length(df.output) >=1)
  {
    
    # Sorting 
    df.output = sort(df.output[!is.na(df.output)], na.last=TRUE, decreasing=FALSE) 
    
    for(i in 1:length(df.output))
    {
      # set new parameter
      if(!is.null(plot.type))
      {
        carpools.raw.genes(untreated.list = list(dataset[[1]], dataset[[2]]), treated.list = list(dataset[[3]],dataset[[4]]), genes=df.output[i], namecolumn=namecolumn, fullmatchcolumn=fullmatchcolumn, do.plot=TRUE, put.names=put.names, type=plot.type)
        
      }
      else
      {
        # Standard plot of Foldchange + Z-Ratio
        carpools.raw.genes(untreated.list = list(dataset[[1]], dataset[[2]]), treated.list = list(dataset[[3]],dataset[[4]]), genes=df.output[i], namecolumn=namecolumn, fullmatchcolumn=fullmatchcolumn, do.plot=TRUE, put.names=put.names, type="foldchange")
        #carpools.raw.genes(untreated.list = list(dataset[[1]], dataset[[2]]), treated.list = list(dataset[[3]],dataset[[4]]), genes=df.output[i], namecolumn=namecolumn, fullmatchcolumn=fullmatchcolumn, do.plot=TRUE, put.names=put.names, type="z-score", sgrna.data=sgrna.data)
        #carpools.raw.genes(untreated.list = list(dataset[[1]], dataset[[2]]), treated.list = list(dataset[[3]],dataset[[4]]), genes=df.output[i], namecolumn=namecolumn, fullmatchcolumn=fullmatchcolumn, do.plot=TRUE, put.names=put.names, type="z-ratio", sgrna.data=sgrna.data)
        carpools.raw.genes(untreated.list = list(dataset[[1]], dataset[[2]]), treated.list = list(dataset[[3]],dataset[[4]]), genes=df.output[i], namecolumn=namecolumn, fullmatchcolumn=fullmatchcolumn, do.plot=TRUE, put.names=put.names, type="foldchange.vioplot", controls.target=controls.target, controls.nontarget=controls.nontarget)
        
      }
    }
    
    return(df.output)
  }
  par(par.old)
    ###### End of creating automatic plots for all hits if labelgenes == NULL

  
  
  
}