## package: megaptera
## author: Christoph Heibl (at gmx.net)
## last update: 2014-04-24

check.Clades <- function(phy, tax){
  
  tax <- tax[, names(tax) != "synonyms"]
  tax <- tax[tax$spec %in% phy$tip.label,]
  
  mp <- function(rank, phy, tax){
    obj <- NULL
    for ( i in unique(tax[, rank]) ){
      taxSet <- tax$spec[tax[, rank] == i]
      s <- noi(phy, taxSet)
      phySet <- descendants(phy, s, labels = TRUE)
      id <- phySet %in% taxSet
      if ( all(id) ) {
        mp <- TRUE
        frac <- 0
      } else {
        mp <- FALSE
        frac <- length(id[!id])/length(id)
      }
      obj <- rbind(obj, c(i, s, mp, frac))
    }
    cbind(rank, obj)
  }
  ranks <- names(tax)[names(tax) != "spec"]
  obj <- lapply(ranks, mp, phy = phy, tax = tax)
  obj <- do.call(rbind, obj)
  obj <- data.frame(obj, stringsAsFactors = FALSE)
  names(obj) <- c("rank", "taxon", "node", "mp", "prop")
  obj$mp <- as.logical(obj$mp)
  obj$prop <- as.numeric(obj$prop)
  obj  
}