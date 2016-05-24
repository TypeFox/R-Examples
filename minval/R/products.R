# products | Identifies the reactants for a stoichometric reaction
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

products <- function(reaction){
  # Identifies if stoichiometric reaction is not reversible. In this case:
  if (grepl("<=>",reaction)){
    products <- unlist(strsplit(reaction,"[[:blank:]]*<=>[[:blank:]]*"))
  }else {
    products <- unlist(strsplit(reaction,"[[:blank:]]*=>[[:blank:]]*"))[2]
  }

  # Split independient reactants
  products <- unlist(strsplit(products,"[[:blank:]]\\+[[:blank:]]"))
  products <- gsub("^[[:blank:]]*","",products)
  products <- gsub("[[:blank:]]*$","",products)
  # Use a regex to extract stoichiometric coefficients and separate the metabolite name
  products <- gsub("^[[:digit:]][[:graph:]]*[[:blank:]]+","",products)
  
  return(products[!is.na(products)])
}

