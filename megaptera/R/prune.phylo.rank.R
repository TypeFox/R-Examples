# LAST UPDATE: 2013-11-21

prune.phylo.rank <- function(phy, tax, rank = "gen", add.nb = TRUE, 
                             quiet = FALSE){
  
  ## Expand taxonomy table if column 'species' is missing.
  ## This happens when GenBank was searched for genera.
  ## --------------------------------------------------
  if ( is.null(tax$spec) ){
    spec <- phy$tip.label
    gen <- levels(tax$genus)[tax$genus]
    new.tax <- data.frame()
    for ( i in seq_along(gen) ){
      id <- grep(paste(gen[i], "_", sep = ""), spec)
      if ( length(id) > 0 ){
        new.tax <- rbind(new.tax, 
                         data.frame(spec = spec[id], tax[i, ], 
                                    row.names = NULL))
      }
    } # end of FOR-loop
    tax <- new.tax
  }
  
  ## taxonomy
  ## --------
  tax <- tax[which(tax$spec %in% phy$tip.label), ]
  
  rank <- split(tax$spec, tax[rank])
  
  for ( i in seq_along(rank) ){
    print(i)
    id <- which(phy$tip.label %in% rank[[i]])
    cn <- names(rank)[i]
    if ( length(id) == 1 ){
      phy$tip.label[id] <- ifelse(add.nb,  
                                  paste(cn, "-", length(id), "specs."), cn)
    }
    else {
      if ( is.monophyletic(phy, id) ) {
        phy$tip.label[id] <- ifelse(add.nb,  paste(cn, "-", length(id), "specs."), cn)
        phy <- drop.tip(phy, id[-1])
        # TO DO: set edge length to crown group
      }
      else {
        if ( !quiet ) message("to bad: ", cn, " is not monophyletic")
      } 
    }
  }
  return(phy)
}