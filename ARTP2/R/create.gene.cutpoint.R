
create.gene.cutpoint <- function(pathway, rs, options){
  
  pd <- pathway[pathway$SNP %in% rs, ]
  GeneInGroup <- unique(pd$Gene)
  ngene <- length(GeneInGroup)
  GeneCutPoint <- list()
  vGeneCutPoint <- NULL
  GeneCutPointStartEnd <- matrix(NA, ngene, 2)
  rownames(GeneCutPointStartEnd) <- GeneInGroup
  colnames(GeneCutPointStartEnd) <- c("Start", "End")
  inspect.snp.n <- options$inspect.snp.n
  inspect.snp.percent <- options$inspect.snp.percent
  
  for(g in 1:ngene){
    gene <- GeneInGroup[g]
    gene.size <- sum(pd$Gene == gene)
    if(inspect.snp.percent == 0){
      cp <- 1:min(inspect.snp.n, gene.size)
    }else{
      step <- floor(gene.size * inspect.snp.percent)
      step <- max(1, step)
      up.limit <- min(gene.size, inspect.snp.n * step)
      up.limit <- max(up.limit, step)
      cp <- seq(step, up.limit, by = step)
    }
    GeneCutPoint[[g]] <- cp
    GeneCutPointStartEnd[gene, "Start"] <- length(vGeneCutPoint) + 1
    vGeneCutPoint <- c(vGeneCutPoint, cp)
    GeneCutPointStartEnd[gene, "End"] <- length(vGeneCutPoint)
  }
  
  list(vGeneCutPoint = vGeneCutPoint, 
       GeneCutPointStartEnd = GeneCutPointStartEnd)
  
}


