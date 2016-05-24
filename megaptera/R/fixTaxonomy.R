# PACKAGE: megaptera
# CALLED BY: user
# AUTHOR: Christoph Heibl
# LAST UPDATE: 2014-11-05

fixTaxonomy <- function(tax, auto = FALSE, ignore = "synonyms"){
  
  ## enforce SQL-compatible attribute names
  tax <- sqlTaxonomyHeader(tax)
  
  ## drop attribute synonym
  tax$synonyms <- NULL
  
  ## adjust order: rightmost column must be the highest rank
  ## -------------------------------------------------------
  if ( which(names(tax) %in% "spec") > which(names(tax) %in% "gen") )
    tax <- tax[, ncol(tax):1]
  
  ## sometimes no genus name is returned by NCBI
  ## -------------------------------------------
  id <- tax$gen == "-"
  if ( any(id) ){
    gen <- strip.spec(tax$spec[id])
    levels(tax$gen) <- union(levels(tax$gen), gen)
    tax$gen[id] <- gen
  }
  
  for ( i in (ncol(tax) - 1):1 ){
#     for ( i in 15:4 ){
    cat("Taxonomic rank:", names(tax)[i])
    xx <- unique(tax[, ncol(tax):i])
    d <- xx[, ncol(xx)]
    d <- d[!d %in% c("-", "incertae sedis")]
    dd <- duplicated(d)
    if ( any(dd) ){
      ## identify and loop over duplicate taxa
      d <- unique(d[dd])
      for ( j in d){
        cat("\nProblem found in: rank '", names(tax)[i], "', taxon '", j, "'", sep = "") 
        ## identify ambiguous higher rank and loop over them
        y <- xx[xx[, ncol(xx)] == j, ]
        y <- apply(y, 2, unique)
        y <- y[sapply(y, length) > 1]
        cat("\n", length(y), "higher ranks are ambiguous")
        for ( k in names(y) ){
          cat("\nRank '", k ,"': please choose one of ...", sep = "")
          z <- y[[k]]
          cat(paste("\n (", 1:length(z), ") '", z, "'", sep = ""))
          if ( auto ){
            ind <- sample(seq_along(z), 1) ## random
          } else {
            ind <- as.numeric(readline()) ## interactive
          }
          cat(" your choice: '", z[ind], "'\n", sep = "")
          tax[tax[, i] == j, k] <- z[ind]
        } # end of FOR-loop over k
      } # end of FOR-loop over j
    } else {
      cat(" ... OK\n")
      d <- NULL
    }
  } # end of FOR-loop over i
  unique(tax[, ncol(tax):1])
}

