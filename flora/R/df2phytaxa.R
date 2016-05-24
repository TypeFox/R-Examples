#' Phylomatic format
#'
#' Convert the results of get.taxa() to the phylomatic sample format
#' 
#' @param taxa A data frame with columns named family, genus, and species.
#' @param uppercase logical. Should the function capitalize first letters?
#' @export
df2phytaxa <- function(taxa, uppercase = TRUE) {
  if (!all(c("family", "search.str") %in% colnames(taxa)) || !inherits(taxa, "data.frame")) {
    stop("taxa is not in the correct format. It should be a data regurgitated by flora::get.taxa")
  }
  res <- vector()
  for (i in seq_len(nrow(taxa))) {
    if (taxa[i, "notes"] == "not found") next
    family <- taxa[i, "family"]
    search.str <- taxa[i, "search.str"]
    taxon.rank <- taxa[i, "taxon.rank"]
    ident <- !grepl("sp\\.", search.str)
    if (taxon.rank == "family" & ident) {
      taxon <- family
    }
    if (taxon.rank == "genus" | (taxon.rank == "family" & !ident)) {
      taxon <- paste(family, gsub(" ", "_", search.str), sep = "/")
    }
    if (taxon.rank == "species" | (taxon.rank == "genus" & !ident)) {
      spp <- taxa[i, "search.str"]
      genus <- strsplit(spp, " ")[[1]][1]
      taxon <- paste(family, genus, gsub(" ", "_", spp), sep = "/")
    }
    res[length(res) + 1] <- taxon
  }
  if (uppercase) {
    res
  } else {
    tolower(res)
  }
}