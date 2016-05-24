# compartments
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

compartments <- function(metabolites){
  # Compartment main function
  compartments <- unique(unlist(regmatches(metabolites, gregexpr("\\[[[:alnum:]]*(\\_)?[[:alnum:]]*\\]$", metabolites))))
  compartments <- gsub("\\[","",compartments)
  compartments <- gsub("\\]","",compartments)
  return(compartments)
}