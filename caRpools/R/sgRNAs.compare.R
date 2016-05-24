sgRNAs.compare=function(deseq.file, mageck.file, genes=NULL, namecolumn=1, fullmatchcolumn=2, extractpattern=expression("^(.+?)_.+"), pval=0.05  )
{
  deseq.file = as.data.frame(deseq.file$sgRNA)
  mageck.file.genes = mageck.file$genes
  mageck.file = as.data.frame(mageck.file$sgRNA) #load.file(paste(mageck.file, "MAGeCK_sgrna_summary.tab", sep="_") )

 
  # Aim:
  # Load all sgRNA files and create two data frames
  # first dataframe has genes + number of significant sgRNAs for all three methods
  # second df has sgRNA + whether they were significant or not
  
  # both can be used for enrich plotting e.g. plot.raw.genes
  
  if(is.null(genes))
  {
    compare.genes = deseq.file$genes
    compare.sgRNA = rownames(deseq.file)
  } else if (!is.null(genes))
  {
    compare.genes = genes
    compare.sgRNA = rownames(deseq.file[deseq.file$genes == genes,])
  }
  
  # span new data.frame
  # use DESeq2 data for creating data.frame
  data.return = data.frame(
    sgRNA = compare.sgRNA,
    genes = compare.genes,
    deseq.padj = 1,
    deseq.sig = "no",
    mageck.fdr.pos = 1,
    mageck.sig.pos = "no",
    mageck.fdr.neg = 1,
    mageck.sig.neg = "no",
    stringsAsFactors = FALSE)

  # structure DESEQ2
  # rownames baseMean  log2FoldChange	lfcSE	stat	pvalue	padj	genes  
  # remove NAs by setting them to 1
  deseq.file$padj = apply(deseq.file, 1, function(x) {
    if(is.finite(as.numeric(x["padj"])))
    {return(as.numeric(x["padj"]))}
    else
    {return(1)}
  })


  
  for(i in 1:nrow(deseq.file))
  {
    if(as.numeric(deseq.file[i,"padj"]) <= pval)
    {
      data.return[data.return$sgRNA == deseq.file[i,"sgRNA"], "deseq.padj"] = deseq.file[i,"padj"]
      
      if(deseq.file[i,"log2FoldChange"]>0)
      {
        data.return[data.return$sgRNA == deseq.file[i,"sgRNA"], "deseq.sig"] = "enriched"
      }
        else if(deseq.file[i,"log2FoldChange"]<0)
        {
          data.return[data.return$sgRNA == deseq.file[i,"sgRNA"], "deseq.sig"] = "depleted"
        }
    }
    
  }

  #print("DESeq2 finished")

  # structure Mageck
  # sgrna  Gene	control_count	treatment_count	control_mean	treat_mean	control_var	adj_var	score	p.low	p.high	p.twosided	FDR	high_in_treatment
  
#   mageck.file$p.low = apply(mageck.file, 1, function(x) {
#     if(is.finite(as.numeric(x["p.low"])))
#     {return(as.numeric(x["p.low"]))}
#     else
#     {return(1)}
#   })
#   mageck.file$p.high = apply(mageck.file, 1, function(x) {
#     if(is.finite(as.numeric(x["p.high"])))
#     {return(as.numeric(x["p.high"]))}
#     else
#     {return(1)}
#   })
  
#   for(i in 1:nrow(mageck.file))
#   {
#     #enriched
#     if(as.numeric(mageck.file[i,"p.high"]) <= pval)
#     {
#       data.return[data.return$sgRNA == mageck.file[i,"sgrna"], "mageck.fdr.pos"] = mageck.file[i,"p.high"]
#       data.return[data.return$sgRNA == mageck.file[i,"sgrna"], "mageck.sig.pos"] = "enriched"
#       data.return[data.return$sgRNA == mageck.file[i,"sgrna"], "mageck.fdr.neg"] = mageck.file[i,"p.low"]
#     }
#     
#     #depleted
#     if(as.numeric(mageck.file[i,"p.low"]) <= pval)
#     {
#       data.return[data.return$sgRNA == mageck.file[i,"sgrna"], "mageck.fdr.neg"] = mageck.file[i,"p.low"]
#       data.return[data.return$sgRNA == mageck.file[i,"sgrna"], "mageck.sig.neg"] = "depleted"
#       data.return[data.return$sgRNA == mageck.file[i,"sgrna"], "mageck.fdr.pos"] = mageck.file[i,"p.high"]  
#     }
#     
#   }
#   
  #print("MAGeCK finished")

  ## get information for certain genes if necessary

  # return for aggregated gene data table
  sgRNA.number.deseq = aggregate.data.frame(data.return$deseq.padj,list(data.return$genes),function(x) length(x[as.numeric(x) != 1]))
  #sgRNA.number.mageck = aggregate.data.frame(data.return$mageck.fdr.pos,list(data.return$genes),function(x) length(x[as.numeric(x) != 1]))
  #sgRNA.number.mageck2 = aggregate.data.frame(data.return$mageck.fdr.neg,list(data.return$genes),function(x) length(x[as.numeric(x) != 1]))
  
  #sgRNA.number.mageck$x = sgRNA.number.mageck$x + sgRNA.number.mageck2$x
  

  data.return2 = data.frame(
    genes = data.return$genes,
    stringsAsFactors=FALSE)
  
  data.return2 = aggregate.data.frame(data.return2, list(data.return2$genes), function(x) return(NA))
  
  data.return2$x=NULL
  data.return2$genes = data.return2$Group.1
  data.return2$Group.1 = NULL


data.return2$deseg.sgRNA = apply(data.return2, 1, function(x) {return(sgRNA.number.deseq[sgRNA.number.deseq$Group.1 == as.character(x["genes"]), "x"])} )
 
#data.return2$mageck.sgRNA = apply(data.return2, 1, function(x) {return(sgRNA.number.mageck[sgRNA.number.mageck$Group.1 == as.character(x["genes"]), "x"])} )
  
data.return2$mageck.sgRNA.pos.good = apply(data.return2, 1, function(x) {return(mageck.file.genes[mageck.file.genes[,"genes"] == as.character(x["genes"]), "sgrna.pos.good"])} )
data.return2$mageck.sgRNA.neg.good = apply(data.return2, 1, function(x) {return(mageck.file.genes[mageck.file.genes[,"genes"] == as.character(x["genes"]), "sgrna.neg.good"])} )


# return as whole dataset that can be used for plot.raw.genes
  return(list(dataset.table = data.return, compare.table = data.return2))
  
  
}

