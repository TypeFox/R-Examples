# metabolites
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

metabolites <- function(reactionList , woCompartment = FALSE){
  # Extract and return the unique metabolites list
  mets <- function(reaction){
    if (grepl("<=>",reaction)){
      metabolites <- unlist(strsplit(reaction,"<=>"))
    } else {
      metabolites <- unlist(strsplit(reaction,"=>"))
    }
    metabolites <- unlist(strsplit(metabolites,"[[:blank:]]\\+[[:blank:]]"))
    # Use a regex to extract stoichiometric coefficients and separate the metabolite name
    metabolites <- gsub("^[[:blank:]]*","",metabolites)
    metabolites <- gsub("[[:blank:]]*$","",metabolites)
    metabolites <- gsub("^[[:digit:]][[:graph:]]*[[:blank:]]","",metabolites)
  }
  metabolites <- as.vector(unique(unlist(sapply(reactionList, mets))))
  if (woCompartment == TRUE){
    return(unique(.metname(metabolites)))
  } else{
    return(metabolites)
  }
}
