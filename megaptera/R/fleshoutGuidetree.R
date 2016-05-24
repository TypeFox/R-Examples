fleshoutGuidetree <- function(phy, tax){
  
  if ( !inherits(phy, "phylo") )
    stop("'phy' is not of class 'phylo'")
  
  ## data frame containing species missing from phy
  ## ----------------------------------------------
  add <- tax[!tax$spec %in% phy$tip.label, ]
  cat("\nNumber of species to add:", nrow(add))
  gen.phy <- strip.spec(phy$tip.label)
  
  ## add via genus
  #for ( i in add$spec[add$gen %in% gen.phy]){
  for ( i in union(add$spec, NULL) ){ # union returns character!
    cat("\n.. ", i)
    phy <- addTip(phy, tip = i, tax = tax)
  }
  phy
}