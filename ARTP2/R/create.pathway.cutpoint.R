
create.pathway.cutpoint <- function(pathway, options){
  
  ngene <- length(unique(pathway$Gene))
  inspect.gene.n <- options$inspect.gene.n
  inspect.gene.percent <- options$inspect.gene.percent
  
  step <- floor(ngene * inspect.gene.percent)
  step <- max(1, step)
  up.limit <- min(ngene, inspect.gene.n * step)
  up.limit <- max(up.limit, step)
  seq(step, up.limit, by = step)
  
}
