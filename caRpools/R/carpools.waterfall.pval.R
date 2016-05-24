carpools.waterfall.pval = function(type=NULL,dataset=NULL, pval=0.05, mageck.type="pos", log=TRUE)
{
  # type 
  # can be deseq2, mageck, wilcox
  
  # dataset
  # must be analysis dataset from stat.deseq2, stat.wilcox, stat.mageck
  pval.old=pval
  
  if(is.null(type))
  {stop("No analysis type selected. Please use either type='deseq2', type='mageck' or type='wilcox'.")}
  else if(type=="wilcox")
  {
    # use wilcoxon information
    plot.data = dataset$p.value
    
  } else if (type=="deseq2")
  {
    # use DESeq2 information
    plot.data = dataset$genes[,"padj"]
    
  } else if(type=="mageck")
  {
    # use MAGeCK information
    plot.data = dataset$genes[,mageck.type]
  }
  
  # transform data to log10 scale
  if(identical(log,TRUE))
  {
    
    # transform to log10 scale
    plot.data = sapply(plot.data, function(x)
      {
        if(is.finite(-log10(as.numeric(x))))
        {
          return(-log10(as.numeric(x)))}
        else return(1)
    })
    
    pval = -log10(pval)
    
    # Sort according to pvalue
    plot.data = sort(plot.data, decreasing=TRUE)
    
    # make plot rdy
    plot.col = sapply(plot.data, function(x)
    {
      if(as.numeric(x) >= pval)
      {return("red")}
      else { return("black")}
    })
    ylab = "-log10 adjusted p-value"
  }
  else
  {
    # Sort according to pvalue
    plot.data = sort(plot.data, decreasing=FALSE)
    
    # make plot rdy
    plot.col = sapply(plot.data, function(x)
    {
      if(as.numeric(x) <= pval)
      {return("red")}
      else { return("black")}
    })
    ylab = "adjusted p-value"
  }
  
  
  
  plot(1:length(plot.data), plot.data, main="p-value Distribution", ylab = ylab, xlab="Gene",
       col= plot.col,
       pch = 16
       )
  # plot pval line
  lines(1:length(plot.data), rep(pval, length(plot.data)), col="black", lty="dashed", lwd=2)
  legend("topright", legend=paste("p-value:",pval.old, sep=" "))
  
  
}