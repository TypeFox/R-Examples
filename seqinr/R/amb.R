amb <- function(base, forceToLower = TRUE, checkBase = TRUE,
IUPAC = s2c("acgturymkswbdhvn"), u2t = TRUE){
  if(missing(base)) return(IUPAC)
  if(!is.character(base)) stop("Character expected")
  if(nchar(base) != 1) stop("Single character expected")
  if(forceToLower) base <- tolower(base)
  if(checkBase) if(!(base %in% IUPAC)) stop("IUPAC base expected")

  if(base == "r") return(s2c("ag")) # puRine
  if(base == "y") return(s2c("ct")) # pYrimidine

  if(base == "m") return(s2c("ac")) # aMino
  if(base == "k") return(s2c("gt")) # Keto

  if(base == "s") return(s2c("cg")) # Strong (3 H bonds)
  if(base == "w") return(s2c("at")) # Weak (2 H bonds)

  if(base == "b") return(s2c("cgt")) # Not A (A->B)
  if(base == "d") return(s2c("agt")) # Not C (C->D)
  if(base == "h") return(s2c("act")) # Not G (G->H)
  if(base == "v") return(s2c("acg")) # Not T (T->V)

  if(base == "n") return(s2c("acgt")) # aNy base
#
# Uracil case: Uracil in RNA instead of Thymine in DNA
#  
  if(base == "u") if(u2t) return("t") else return("u")

  return(base)
}