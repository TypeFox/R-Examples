carpools.reads.genedesigns = function(dataset, namecolumn=1, fullmatchcolumn=2, title="Read Count", xlab="Percentage of sgRNAs present", ylab="Number of Genes", agg.function=sum, extractpattern=expression("^(.+?)_.+"), col = rgb(0, 0, 0, alpha = 0.65))
{
  # plot fullmatchreads per gene on X axis
  # plot number of fullmatchdesigns on Y axis

  
  # get gene names
  rownames(dataset) = dataset[,namecolumn]
  gene.names = sub(extractpattern,"\\1",dataset[,namecolumn],perl=TRUE)
  dataset[,namecolumn] = gene.names
  
  dataset$genes = gene.names
  
  # get readcount
  dataset.genes = aggregate(dataset, list(dataset$genes),
                        FUN=function(x){
                          if(is.numeric(x)){
                            agg.function(x)
                          }else{
                            x[1]
                          }
                        }) 
  
  # get number of max designs per gene
  design.max = aggregate(dataset[,fullmatchcolumn], by=list(dataset$genes), function(x) return(length(x)))

  design.present.all = apply(dataset, 1, function(x) if(as.numeric(x[fullmatchcolumn])>0) {return(as.numeric(1))} else { return(as.numeric(0)) })
  
  design.present = aggregate(design.present.all, by=list(dataset$genes), function(x) return(length(x[x==1])))

  design.present.final = data.frame(gene = design.max$Group.1,
                                    design.all = design.max$x,
                                    design.present = design.present$x,
                                    design.fraction = (design.present$x / design.max$x)*100,
                                    stringsAsFactors=FALSE)
  
  # add read counts
  hist(design.present.final$design.fraction, breaks=50, xlab=xlab, ylab=ylab, main=title, border=FALSE, col=col)
  
  
}

