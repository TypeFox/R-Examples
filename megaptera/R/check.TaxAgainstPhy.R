check.TaxAgainstPhy <- function(phy, tax, rank = NULL){
  
  if ( !inherits(phy, "phylo") )
    stop("argument 'phy' is not of class 'phylo'")
  
  ## intersect tax and phy
  ## ---------------------
  tax <- tax[tax$spec %in% phy$tip.label, ]
  phy <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% tax$spec])
  
  x <- split(tax$spec, tax[rank])
  
  obj <- data.frame(taxon = names(x))
  obj <- cbind(obj, 
               topo = sapply(x, is.monophyletic, phy = phy))
  obj$topo[obj$topo] <- "mono"
  obj <- cbind(obj, size = sapply(x, length))
  obj
}