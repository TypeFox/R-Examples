"aa123" <-
function (aa) {
  if(any(nchar(aa)!=1))
    stop("Provide a character vector of individual 1-letter aminoacid codes")
  
  # convert one-letter IUPAC amino-acid code into
  # three-letter PDB style, for instance "A" into "ALA".
  
  aa1 <- c("-","X",
           "A","C","D","E","F","G",
           "H","I","K","L","M","N","P","Q",
           "R","S","T","V","W","Y")
  aa3 <- c("---","UNK",
           "ALA", "CYS", "ASP", "GLU", "PHE", "GLY",
           "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN",
           "ARG", "SER", "THR", "VAL", "TRP", "TYR")

  convert <- function(x) {
    if(is.na(x)) return(NA)
    if (all(x != aa1)) {
      warning("Unknown one letter code for aminoacid")
      return("UNK")
    }
    else {
      return(aa3[which(x == aa1)])
    }
  }
  return(as.vector(unlist(sapply(aa, convert))))
}

