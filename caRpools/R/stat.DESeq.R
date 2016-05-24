stat.DESeq=function(untreated.list,treated.list,namecolumn=1, fullmatchcolumn=2, agg.function=sum, extractpattern=expression("^(.+?)_.+"), sorting=FALSE, sgRNA.pval = 0.01, filename.deseq="data", fitType="parametric", p.adjust="holm"){

  #requireNamespace(DESeq2)
  non=lapply(
    c(untreated.list,
      treated.list),
    FUN=function(x) stopifnot(
      identical(
        x[,namecolumn],
        treated.list[[1]][,namecolumn]
      )
    )
  )

  # FIRST run for sgRNA enrichment/depletion
  # THEN run for gene effect
  # write table with sgRNA infos that will be loaded later for in-depth analysis
  # add $sgRNAs to res for number of significant sgRNAs
  
  
  designs=treated.list[[1]][,namecolumn]    
  
  countdata=do.call(
    "cbind.data.frame",
    lapply(
      c(
        untreated.list,
        treated.list
      ),
      FUN=function(x){
        as.numeric(x[,fullmatchcolumn])
      }
    )
  )
  
  # Run for testing sgRNA effects
  row.names(countdata)=treated.list[[1]][,namecolumn]

  colnames(countdata)=c(rep("CTRL",length(untreated.list)),rep("TREAT",length(treated.list)))

  coldata=cbind.data.frame(condition = c(rep("control",length(untreated.list)),rep("treated",length(treated.list))))
  row.names(coldata)=c(names(untreated.list),names(treated.list))

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = countdata,
                                     colData = coldata,
                                design =  ~ condition)
  
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- DESeq2::estimateDispersions(dds, fitType="parametric")
  des.sgRNA <- DESeq2::nbinomWaldTest(dds)
  
  #plotDispEsts(dds)

  res.sgRNA <- DESeq2::results(des.sgRNA)
  
  # now run for gene enrichment
  # aggregate data first
  
  gene.names = sub(extractpattern,"\\1",treated.list[[1]][,namecolumn],perl=TRUE)  
  gene.names.agg=aggregate(gene.names,list(gene.names),function(x) x[1])$Group.1
  
  countdata=do.call(
    "cbind.data.frame",
    lapply(
      c(
        untreated.list,
        treated.list
      ),
      FUN=function(x){as.integer(
        aggregate(
          x[,fullmatchcolumn],
          list(gene.names),
          agg.function
        )$x
      )}
    )
  )
  
  row.names(countdata)=gene.names.agg
  
  colnames(countdata)=c(rep("CTRL",length(untreated.list)),rep("TREAT",length(treated.list)))
  
  coldata=cbind.data.frame(condition = c(rep("control",length(untreated.list)),rep("treated",length(treated.list))))
  row.names(coldata)=c(names(untreated.list),names(treated.list))
  dds.genes <- DESeq2::DESeqDataSetFromMatrix(countData = countdata,
                                colData = coldata,
                                design =  ~ condition)
    
  # TRY CLEAN UP FOR OUTLIERS
  dds.genes <- DESeq2::estimateSizeFactors(dds.genes)
  dds.genes <- DESeq2::estimateDispersions(dds.genes, fitType=fitType)
  des.genes <- DESeq2::nbinomWaldTest(dds.genes)
  
  #plotDispEsts(des.genes)
  res.genes <- DESeq2::results(des.genes, pAdjustMethod = p.adjust)
  
  # Combine
  res.sgRNA$genes = sub(extractpattern,"\\1",row.names(res.sgRNA),perl=TRUE)
  res.sgRNA$sgRNA = row.names(res.sgRNA)
  
  # count number of sgRNAs for a given gene
  sgRNA.number = aggregate.data.frame(res.sgRNA$padj,list(res.sgRNA$genes),function(x) length(x[as.character(x) <= sgRNA.pval]))
  
  res.genes$genes = row.names(res.genes)
  res.genes$sgRNA = apply(as.data.frame(res.genes), 1, function(x) {
                   return(sgRNA.number[sgRNA.number[,1] == x["genes"],2])
  })
  
  
  
  # check for NAs in padjusted and pvalue
  res.genes[is.na(res.genes[,"pvalue"]), "pvalue"] <- 1
  res.genes[is.na(res.genes[,"padj"]), "padj"] <- 1
  res.genes[is.na(res.genes[,"baseMean"]), "baseMean"] <- 0
  res.genes[is.na(res.genes[,"log2FoldChange"]), "log2FoldChange"] <- 0
  res.genes[is.na(res.genes[,"stat"]), "stat"] <- 0
  ### sgRNAs
  res.sgRNA[is.na(res.sgRNA[,"pvalue"]), "pvalue"] <- 1
  res.sgRNA[is.na(res.sgRNA[,"padj"]), "padj"] <- 1
  res.sgRNA[is.na(res.sgRNA[,"baseMean"]), "baseMean"] <- 0
  res.sgRNA[is.na(res.sgRNA[,"log2FoldChange"]), "log2FoldChange"] <- 0
  res.sgRNA[is.na(res.sgRNA[,"stat"]), "stat"] <- 0
  
  # Write sgRNA table
  write.table(res.sgRNA, file=paste(filename.deseq, "DESeq2_sgRNA.tab", sep="_"), col.names=TRUE, quote=FALSE, sep="\t", row.names=TRUE)
  
  if(sorting==TRUE){
    res.genes = res.genes[order(res.genes$padj, na.last=TRUE),]
    res.sgRNA = res.sgRNA[order(res.sgRNA$padj, na.last=TRUE),]
  }
  
# Return list
return(list(genes = res.genes, sgRNA = res.sgRNA, deseq.genes = des.genes, deseq.sgRNA = des.sgRNA))

  
}
