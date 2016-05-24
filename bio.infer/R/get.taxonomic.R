get.taxonomic <- function(bcnt) {

  # Global variables for java version: parse.list, fulltab, dupe.list, bcnt
  # Presence of non-null tname.new, and isel determine entry points

  names0 <- names(bcnt)
  get.tax.env <- environment(NULL)
  tlevs <- load.itis(get.tax.env)

  # check for duplicate name with ITIS
  for (i in 1:length(names0)) {
    if (names0[i] %in% tlevs) {
      names0[i] <- paste(names0[i], "orig", sep = ".")
    }
  }
  names(bcnt) <- names0
  
  tname <- get.taxon.names(bcnt)
  df.parse <- parse.taxon.name(tname)
  
  parse.list <- get.valid.names(df.parse, get.tax.env)

  if (length(parse.list) > 1) {
    if (nrow(parse.list[[2]]) > 0) {
      
      # resolve records with multiple valid ITIS entries
      parse.list <- resolve.mult(parse.list, get.tax.env)

      if (nrow(parse.list[[2]]) > 0) {
    # prompt user to correct misspellings
        tname.new <- correct.taxanames(sort(unique(parse.list[[2]][,2])),
                                     get.tax.env)

    # incorporate corrected spellings
        parse.list <- incorp.correct(tname.new, parse.list)
      }
    }
  }

  fulltab <- make.fulltab1(parse.list[[1]], get.tax.env)

  # identify and fix duplicate entries
  dupe.list <- locate.dupes(fulltab)
  if (length(dupe.list$isav) > 0) {
    dupe.sel <- get.dupe.sel(dupe.list$sumstr)
    fulltab <- remove.dupes(fulltab, dupe.list, dupe.sel)
  }

  # assign species names
  finaltab <- make.species(parse.list[[1]], fulltab)
  output.tax.table(finaltab, tlevs)

  bcnt.new <- merge(bcnt, finaltab, by.x = names0[2], by.y = "taxaname.orig")
  return(bcnt.new[, c(names0,tlevs, "SPECIES")])
}
