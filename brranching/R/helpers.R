phylo_base <- "http://phylodiversity.net/phylomatic/pmws"

collapse_double_root <- function(y) {
  temp <- strsplit(y, ")")[[1]]
  double <- c(length(temp) - 1, length(temp))
  tempsplit <- temp[double]
  tempsplit_1 <- strsplit(tempsplit[1], ":")[[1]][2]
  tempsplit_2 <- strsplit(tempsplit[2], ":")[[1]]
  rootlength <- as.numeric(tempsplit_1) +
    as.numeric(strsplit(tempsplit_2[2], ";")[[1]][1])
  newx <- paste(")", tempsplit_2[1], ":", rootlength, ";", sep = "")
  newpre <- gsub("[(]", "", temp[1])
  allelse <- temp[-1]
  allelse <- allelse[setdiff(1:length(allelse), double - 1)]
  allelse <- paste(")", allelse, sep = "")
  paste(newpre, paste(allelse, collapse = ""), newx, sep = "")
}

colldouble <- function(z) {
  if ( class( try( read.tree(text = z), silent = TRUE ) ) %in% 'try-error' ) {
    treephylo <- collapse_double_root(z)
  } else {
    treephylo <- z
  }
  return(treephylo)
}

getnewick <- function(x) {
  tree <- gsub("\n", "", x[[1]])
  read.tree(text = colldouble(tree))
}
