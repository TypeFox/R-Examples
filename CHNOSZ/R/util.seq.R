# CHNOSZ/util.seq.R
# functions to work with sequence data

aminoacids <- function(nchar=1, which=NULL) {
  # return the abbreviations or names of the amino acids
  # the following are all in the same order as thermo$protein
  # the single-letter codes
  aa1 <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
           "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  # the 3-letter codes
  aa3 <- c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
            "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")
  # names of the neutral amino acids
  aaneutral <- c("alanine", "cysteine", "aspartic acid", "glutamic acid", "phenylalanine", 
    "glycine", "histidine", "isoleucine", "lysine", "leucine", "methionine", "asparagine", 
    "proline", "glutamine", "arginine", "serine", "threonine", "valine", "tryptophan", "tyrosine")
  # names of the amino acids with ionized side chains
  aacharged <- c("alanine", "cysteinate", "aspartate", "glutamate", "phenylalanine", 
    "glycine", "histidinium", "isoleucine", "lysinium", "leucine", "methionine", "asparagine", 
    "proline", "glutamine", "argininium", "serine", "threonine", "valine", "tryptophan", "tyrosinate")
  # defaults are in the same order as in thermo$protein
  if(is.null(which)) which <- aa1
  # figure out which amino acids are wanted (can use 1- or 3-letter codes, or neutral names)
  if(all(nchar(which)==1)) iaa <- match(which, aa1)
  else if(all(nchar(which)==3)) iaa <- match(which, aa3)
  else iaa <- match(which, aaneutral)
  # return the desired abbreviations or names
  if(nchar==1) return(aa1[iaa])
  else if(nchar==3) return(aa3[iaa])
  else if(nchar=="") return(aaneutral[iaa])
  else if(nchar=="Z") return(aacharged[iaa])
}

nucleic.formula <- function(nucleic=NULL) {
  # compute the formula, e.g.
  # DNA <- count.aa(list("AGCT", "TTTT"), type="DNA")  # a dataframe of counts
  # nf <- nucleic.formula(DNA)  # a series of formulas
  # !!! this only adds the formulas of the nucleobases; dehydration and phosphorylation are not yet accounted for !!!
  # 20090926 jmd
  letts <- c("A", "C", "G", "T", "U")
  names <- c("adenine", "cytosine", "guanine", "thymine", "uracil")
  # the locations of the letters in the data frame
  i.lett <- match(letts, colnames(nucleic))
  # we'll normally have at least one NA (U or A for DNA or RNA)
  ina <- is.na(i.lett)
  # the species indices of the bases, in the order appearing above
  i.base <- suppressMessages(info(names[!ina], check.it=FALSE))
  # the chemical formula of bases
  f.base <- get("thermo")$obigt$formula[i.base]
  # loop over the base counts
  f.out <- character()
  for(i in 1:nrow(nucleic)) {
    # use makeup() with multipliers and sum=TRUE  20120119 jmd
    f <- as.chemical.formula(makeup(f.base, multiplier=as.numeric(nucleic[i, i.lett[!ina]]), sum=TRUE))
    f.out <- c(f.out, f)
  }
  return(f.out)
}

nucleic.complement <- function(nucleic=NULL, type="DNA") {
  # return the nucleobase complement
  # nucleic.complement(nucleic, "DNA")  # DNA complement
  # nucleic.complement(nucleic, "RNA")  # RNA complement
  # the reference sequence, and its DNA and RNA complements
  ref <- c("A", "C", "G", "T", "U")
  DNA <- c("T", "G", "C", "A", "A")
  RNA <- c("U", "G", "C", "A", "A")
  iref <- match(colnames(nucleic), ref)
  i.base <- which(!is.na(iref))
  colnames(nucleic)[i.base] <- get(type)[iref[i.base]]
  # be nice and re-alphabetize the columns
  o.base <- order(colnames(nucleic)[i.base])
  nucleic <- nucleic[, i.base[o.base], drop=FALSE]
  return(nucleic)
}
