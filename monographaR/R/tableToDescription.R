tableToDescription <-
function(data, filename="species_descriptions.txt") {
  if (class(data) != "data.frame") {
    stop("data must be a data.frame")
  }
  message("Assuming the columns are ordered as: Character_in_description, complement, separator, and the species in the remaining columns")
  cat(" ", fill=T)
  data[,1:3] -> model
  data[,4:ncol(data)] -> spp.data
  colnames(spp.data) -> spp
  cat("Species Descriptions", file=filename)
  spaces <-function(x) return(gsub("(?<=[\\s])\\s*|^\\s+$", "", x, perl=TRUE))
  for (s in 1:length(spp)) {
    cat("\n", file=filename, fill=T, append=T)
    cat(sub(".", " ", spp[s], fixed=T), file=filename, fill=T, append=T)
    cat("\n", file=filename, fill=F, append=T)
    for (k in 1:nrow(data)) {
      model[k,] -> m0
      as.matrix(m0) -> m0
      m0[,3] -> sep0
      m0[,1] -> pre0
      m0[,2] -> suf0
      if (pre0 == "") {
        pre0[] <- "xxxwww"
      }
      if (suf0 == "") {
        suf0[] <- "xxxwww"
      }
      as.character(spp.data[k,s]) -> d0
      paste(pre0, " ", d0, " ", suf0, sep0, " ", sep="") -> s0
      sub("xxxwww ", "", s0) -> s0
      sub("xxxwww ", "", s0) -> s0
      sub(" xxxwww", "", s0) -> s0
      spaces(s0) -> s0
      cat(s0, file=filename, fill=F, append=T)
    }
  }
  if (nchar(filename) > 0) {
    cat("Species descriptions were saved in:")
    cat("\n", getwd())
  }
}
