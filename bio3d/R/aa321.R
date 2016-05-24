"aa321" <-
function (aa) {

  # convert three-letters amino-acid code into
  # one-letter IUPAC code, for instance "ALA" into "A".
  # new residues should be added to through the util/make_aatable.R script

  aa1 <- c("-",   ".",  "X",   bio3d::aa.table$aa1)
  aa3 <- c("---", "---","UNK", bio3d::aa.table$aa3)

    convert <- function(x) {
      if(is.na(x)) return(NA)
      if (all(x != aa3)) {
        warning(paste("Unknown 3-letters code for aminoacid:",x))
        return("X") # mask unk
      }
      else {
        return(aa1[which(x == aa3)])
      }
    }
  return(as.vector(unlist(sapply(aa, convert))))
}


".aa321.na" <-
function (aa) {

  # convert three-letters amino-acid code into
  # one-letter IUPAC code, for instance "ALA" into "A".

  aa1 <- c("-",".","X",
           "C",  "G",  "T",  "A", "U", "I",
           "C",  "G",  "T","  A", "U", "I")
  aa3 <- c("---", "---","UNK",
           "DC", "DG", "DT", "DA", "DU", "DI",
            "C",  "G",  "T",  "A",  "U", "I")
  
    convert <- function(x) {
      if(is.na(x)) return(NA)
      if (all(x != aa3)) {
        warning(paste("Unknown 3-letters code for residue:",x))
        return("X") # mask unk
      }
      else {
        return(aa1[which(x == aa3)])
      }
    }
  return(as.vector(unlist(sapply(aa, convert))))
}

