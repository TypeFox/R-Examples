timeSeq.geneset = function(data, group.label, genesets, gene.names, geneset.names = NULL, 
	reads = NULL, exon.length = NULL, gene.length = NULL, exon.level = TRUE, 
	gene.level = FALSE, p.values = FALSE, n_cores = NULL, iterations = NULL, offset = TRUE){
  
    length.set = length(genesets)
    if (is.null(geneset.names)) geneset.names = paste("set", as.character(1 : length.set), sep = "")

    ##  call timeSeq function calculate the true KLRs.
    model = timeSeq(data, group.label, gene.names, reads, exon.length, gene.length, exon.level, gene.level, p.values, n_cores, iterations, offset)
    NPDE.ratio = model$NPDE.ratio
    PDE.ratio = model$PDE.ratio
    unique.gene.names = unique(gene.names)  
	geneset.NPDE.ratio = numeric(length.set)
	geneset.PDE.ratio = numeric(length.set)	
	for(i in 1 : length.set) {
			geneset.NPDE.ratio[i] = mean(NPDE.ratio[unique.gene.names %in% genesets[[i]]])
			geneset.PDE.ratio[i] = mean(PDE.ratio[unique.gene.names %in% genesets[[i]]])
	}
	
	geneset.NPDE.pvalue = NULL
	geneset.PDE.pvalue = NULL	
    if (p.values) {
    	NPDE.pvalue = model$pvalue[, 2]
    	PDE.pvalue = model$pvalue[, 3]
    	for(i in 1 : length.set) {
    			geneset.NPDE.pvalue[i] = mean(NPDE.pvalue[unique.gene.names %in% genesets[[i]]])
				geneset.PDE.pvalue[i] = mean(PDE.pvalue[unique.gene.names %in% genesets[[i]]])
    	}
    }
	
    out = list(NPDE.ratio = geneset.NPDE.ratio,
			   PDE.ratio = geneset.PDE.ratio,
   		       NPDE.pvalue = geneset.NPDE.pvalue, 
			   PDE.pvalue = geneset.PDE.pvalue,
               geneset.names = geneset.names)
    class(out) = "timeSeq.geneset"
    out
}
