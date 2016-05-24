is.aa <- function (seq, aa.strict = FALSE) {

  #one letter codes for amino acids 
  aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
    "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "J", "X")

  #remove ambiguous amino acids
  if (aa.strict)
    aa <- aa[1:20]

  return(toupper(seq) %in% aa)
}