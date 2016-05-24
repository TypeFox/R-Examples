# reactants | Identifies the reactants for a stoichometric reaction
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

reactants <- function(reaction){
  # Identifies if stoichiometric reaction is not reversible. In this case:
  if (grepl("<=>",reaction)){
    reactants <- unlist(strsplit(reaction,"[[:blank:]]*<=>[[:blank:]]*"))
  }
  # In contrary case:
  else {
    reactants <- unlist(strsplit(reaction,"=>"))[1]
  }
  # Split independient reactants
  reactants <- unlist(strsplit(reactants,"[[:blank:]]*\\+[[:blank:]]*"))
  reactants <- gsub("^[[:blank:]]*","",reactants)
  reactants <- gsub("[[:blank:]]*$","",reactants)
  # Use a regex to extract stoichiometric coefficients and separate the metabolite name
  reactants <- gsub("^[[:digit:]][[:graph:]]*[[:blank:]]","",reactants)
  return(reactants[!is.na(reactants)])
}

