# PACKAGE: megaptera
# CALLED BY: megapteraProj, stepF
# AUTHOR: Christoph Heibl
# LAST UPDATE: 2014-07-28

tax2tree <- function(tax, tip.rank = "spec", ignore = "synonyms"){
  
  ## enforce SQL-compatible attribute names
  tax <- sqlTaxonomyHeader(tax)
  
  ## drop 'ignore' column
  if ( ignore %in% names(tax) ) tax[ignore] <- NULL
  
  ## adjust order: rightmost column must be the highest rank
  ## -------------------------------------------------------
  if ( all(c("spec", "gen") %in% names(tax)) ){
    if ( which(names(tax) %in% "spec") > which(names(tax) %in% "gen") )
      tax <- tax[, ncol(tax):1]
  } else {
    if ( !all(c("fam", "gen") %in% names(tax)) ) stop("no columns spec, gen, fam")
    if ( which(names(tax) %in% "gen") > which(names(tax) %in% "fam") )
      tax <- tax[, ncol(tax):1]
  }
  
  ## ajust to 'tip.rank'
  ## -------------------
  tax <- tax[, grep(paste("^", tip.rank, "$", sep = ""), 
                    names(tax)):ncol(tax), drop = FALSE]
  tax <- unique(tax)
  
  if ( ncol(tax) == 1 | nrow(tax) == 1 ) return(NULL)
  
  ## make columns factors
  ## --------------------
  id <- apply(tax, 2, is.factor)
  id <- names(id)[!id]
  for ( i in id ){
    tax[, i] <- factor(tax[, i])
  }
  
  ## replace missing values by hyphens; this means
  ## these ranks do not exist for the taxa concerned
  ## not sure if this is always reasonable
  ## -------------------------------------
  #tax[is.na(tax)] <- "-"
  #stop("missing values not allowed")
  
  ## delete ranks with one single level
  ## ----------------------------------
  id <- apply(tax, 2, unique)
  id <- sapply(id, length)
  tax <- tax[, id > 1, drop = FALSE]
  
  ## only one rank is available
  if ( ncol(tax) == 1 ){
    gt <- paste(levels(tax[, 1]), collapse = ",")
    gt <- paste("(", gt, ");", sep = "")
    gt <- read.tree(text = gt)
  } else { ## more than one rank is available
    if ( any(is.na(tax)) ) 
      stop("no missing values allowed in 'tax'")
    ff <- do.call(list, tax)
#     ff <- apply(tax, 2, unique) # does not work, when all ranks have same number of taxa
    ff <- sapply(ff, function(x) length(unique(x)))
    ff <- names(ff)[ff > 1]
    ff <- ff[ff != "species"]
    gt <- unique(tax[ff])
    cat("\nguide tree created from taxonomy using ranks:",
        paste("\n - ", ff, sep = ""), "\n")
    ff <- paste("~", paste(rev(ff), collapse = "/"), sep = "")
    gt <- as.phylo(formula(ff), data = gt)
  }
  
  ## duplicate tips are caused by faulty higher ranks
  ## e.g. some Osmia in family Megachilidae and some in 
  ## Megacihlidae
  if ( any(duplicated(gt$tip.label)) )
    stop("guide tree must not contain duplicated tip labels")
  
  gt
}