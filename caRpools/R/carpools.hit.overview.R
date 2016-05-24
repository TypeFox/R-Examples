carpools.hit.overview = function(wilcox=NULL, deseq=NULL, mageck=NULL, cutoff.deseq = 0.001, cutoff.wilcox = 0.05, cutoff.mageck = 0.05, cutoff.override=FALSE, cutoff.hits=NULL, plot.genes="overlapping", type="all")
{
  ## The idea is a plot combining all hits in a single plot.
  # Either all overlapping hits are plotted for enriched and depleted dependent on log2FoldChange on Y-Axis and P-value on X-axis or MAGeCK Rank.
  # By default, the FDR corrected P-value of MAGeCK is used for plotting on the Y-Axis
  
  # Enriched hits are RED, depleted are BLUE, all other GREY
  # Moreover, the P-value cutoff of mageck is plotted as a line
  
  # Call function to gene hit list to get OVERLAPPING genes or all individual genes (if plot.genes != overlapping)
  df.output.enriched = generate.hits(wilcox=wilcox, deseq=deseq, mageck=mageck, type="enriched", cutoff.deseq = cutoff.deseq, cutoff.wilcox = cutoff.wilcox, cutoff.mageck = cutoff.mageck, cutoff.override=cutoff.override, cutoff.hits=cutoff.hits, plot.genes=plot.genes)
  df.output.depleted = generate.hits(wilcox=wilcox, deseq=deseq, mageck=mageck, type="depleted", cutoff.deseq = cutoff.deseq, cutoff.wilcox = cutoff.wilcox, cutoff.mageck = cutoff.mageck, cutoff.override=cutoff.override, cutoff.hits=cutoff.hits, plot.genes=plot.genes)
  

  if(type == "all")
  { par(mfrow=c(1,2)) }
  else {par(mfrow=c(1,1))}
  
    # Moreover we need to generate dataset information
    # All information is passed on by individual datasets from the hit analysis
    
    # Use DESeq2 as provider of the log2Foldchange
    if(!is.null(deseq))
    {
      df.deseq=deseq$genes
    }
    if(!is.null(wilcox))
    {
      df.wilcox=data.frame(p.value = wilcox[,"p.value"], log2fc = log2(wilcox[,"foldchange"]), stringsAsFactors=FALSE)
      df.wilcox$genes = rownames(wilcox)
      
    }
    
  #str(df.deseq[1:10,])
    # Use MAGeCK as provider of p-values
  df.mageck = mageck$genes
  
  #str(df.mageck[,])
    # Now we get the information and create a new dataset using all data
    df.plot = data.frame(
      genes = df.mageck$genes,
      pval.enriched = as.numeric(df.mageck$pos),
      pval.depleted = as.numeric(df.mageck$neg),
      stringsAsFactors=FALSE)
  #str(df.plot[1:10,])  
  
    df.plot$log2fc = apply(df.plot,1, function(x){
      
      return(df.deseq[df.deseq$genes == x["genes"],"log2FoldChange"])
      
    })
    df.plot$deseq = apply(df.plot,1, function(x){
      
      return(df.deseq[df.deseq$genes == x["genes"],"padj"])
      
    })
    df.plot$wilcox = apply(df.plot,1, function(x){
      
      return(df.wilcox[df.wilcox$genes == x["genes"],"p.value"])
      
    })
  df.plot$wilcox.log2fc = apply(df.plot,1, function(x){
    
    return(df.wilcox[df.wilcox$genes == x["genes"],"log2fc"])
    
  })
    
    # now we add color
  
  #print(df.plot[1:10,])
  
  if(type=="all" || type=="enriched")
  {
    df.plot$color.enriched = apply(df.plot, 1, function(x){
      
      ret.col = rgb(211,211,211, 255, maxColorValue=255)
      
      if(length(df.output.enriched) >= 1)
      {
        if(x["genes"] %in% df.output.enriched)
        {ret.col = rgb(217,35,35, 255, maxColorValue=255)}
        else
        {
          if(as.numeric(x["deseq"]) < cutoff.deseq && as.numeric(x["log2fc"]) > 0)
          {
            ret.col = rgb(255,165,0, 255, maxColorValue=255)
          }
          if(as.numeric(x["wilcox"]) < cutoff.wilcox && as.numeric(x["wilcox.log2fc"]) > 0)
          {
            ret.col = rgb(255,165,0, 255, maxColorValue=255)
          }
          if(as.numeric(x["pval.enriched"]) < cutoff.mageck)
          {
            ret.col = rgb(255,165,0, 255, maxColorValue=255)
          }
        }
        
        
      }
      else # no overlapping hits!
      {
        
        if(as.numeric(x["deseq"]) < cutoff.deseq && as.numeric(x["log2fc"]) > 0)
        {
          ret.col = rgb(255,165,0, 255, maxColorValue=255)
        }
        if(as.numeric(x["wilcox"]) < cutoff.wilcox && as.numeric(x["wilcox.log2fc"]) > 0)
        {
          ret.col = rgb(255,165,0, 255, maxColorValue=255)
        }
        if(as.numeric(x["pval.enriched"]) < cutoff.mageck)
        {
          ret.col = rgb(255,165,0, 255, maxColorValue=255)
        }
        
      }
      return(ret.col)
    })

    # Enriched
    
    
    plot(abs(log10(df.plot$pval.enriched)),
         df.plot$log2fc,
         col=df.plot$color.enriched,
         pch=16,
         ylab="DESeq2 log2(foldchange)",
         xlab="MAGeCK enriched -log10(p-value)",
         main="Enriched Genes",
         cex=1.1)
    
    # Generate p-value line
    abline(v=abs(log10(cutoff.mageck)))
    
    # Generate Lebel for overlapping hits
    if(length(df.output.enriched) >= 1)
    {
      # Enriched
      enriched.genes=c(which(df.plot$genes %in% df.output.enriched))
      text(abs(log10(df.plot$pval.enriched))[enriched.genes],df.plot$log2fc[enriched.genes],df.plot$genes[enriched.genes],cex=1,pos=2,offset = 1)
      
    }
    # legend
    legend("topleft",bty = "n", cex = 0.9, pt.bg = c(rgb(217,35,35, 255, maxColorValue=255),rgb(255,165,0, 255, maxColorValue=255),rgb(211,211,211, 255, maxColorValue=255) ), col = c(rgb(217,35,35, 255, maxColorValue=255),rgb(255,165,0, 255, maxColorValue=255),rgb(211,211,211, 255, maxColorValue=255) ), pch = c(16,16,16), legend=c("Overlapping Hit", "Non-Overlapping Hit", "Not Enriched") )
  }
 
  if(type=="all" || type=="depleted")
  {
      df.plot$color.depleted = apply(df.plot, 1, function(x){
        
        ret.col = rgb(211,211,211, 255, maxColorValue=255)
        
        if(length(df.output.depleted) >= 1 )
        {
          if(x["genes"] %in% df.output.depleted)
          {
            ret.col = rgb(46,98,166, 255, maxColorValue=255)
          }
          else
          {
            if(as.numeric(x["deseq"]) < cutoff.deseq && as.numeric(x["log2fc"]) < 0)
            {
              ret.col = rgb(255,165,0, 255, maxColorValue=255)
            }
            if(as.numeric(x["wilcox"]) < cutoff.wilcox && as.numeric(x["wilcox.log2fc"]) < 0)
            {
              ret.col = rgb(255,165,0, 255, maxColorValue=255)
            }
            if(as.numeric(x["pval.depleted"]) < cutoff.mageck)
            {
              ret.col = rgb(255,165,0, 255, maxColorValue=255)
            }
          }
          
        }
        else
        {
          if(as.numeric(x["deseq"]) < cutoff.deseq && as.numeric(x["log2fc"]) < 0)
          {
            ret.col = rgb(255,165,0, 255, maxColorValue=255)
          }
          if(as.numeric(x["wilcox"]) < cutoff.wilcox && as.numeric(x["wilcox.log2fc"]) < 0)
          {
            ret.col = rgb(255,165,0, 255, maxColorValue=255)
          }
          if(as.numeric(x["pval.depleted"]) < cutoff.mageck)
          {
            ret.col = rgb(255,165,0, 255, maxColorValue=255)
          }
          
        }
        return(ret.col)
      })
      
      ## Depleted plot
      plot(abs(log10(df.plot$pval.depleted)),
           df.plot$log2fc,
           col=df.plot$color.depleted,
           pch=16,
           ylab="DESeq2 log2(foldchange)",
           xlab="MAGeCK depleted -log10(p-value)",
           main="Depleted Genes",
           cex=1.1)
      
      # Generate p-value line
      abline(v=abs(log10(cutoff.mageck)))
      
      # Generate Lebel for overlapping hits
    if(length(df.output.depleted) >= 1 )
    {
      # depleted
      depleted.genes=c(which(df.plot$genes %in% df.output.depleted))
      text(abs(log10(df.plot$pval.depleted))[depleted.genes],df.plot$log2fc[depleted.genes],df.plot$genes[depleted.genes],cex=1,pos=2,offset = 1)
      
    }
    # legend
    legend("topleft",bty = "n", cex = 0.9, pt.bg = c(rgb(46,98,166, 255, maxColorValue=255),rgb(255,165,0, 255, maxColorValue=255),rgb(211,211,211, 255, maxColorValue=255) ), col = c(rgb(46,98,166, 255, maxColorValue=255),rgb(255,165,0, 255, maxColorValue=255),rgb(211,211,211, 255, maxColorValue=255) ), pch = c(16,16,16), legend=c("Overlapping Hit", "Non-Overlapping Hit", "Not Depleted") )
    
  }
  
  # set par back to normal
  par(mfrow=c(1,1))
  
}