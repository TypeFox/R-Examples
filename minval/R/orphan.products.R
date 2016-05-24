# orphan.products
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

orphan.products <- function(reactionList, byCompartment=FALSE){
  # Extract all reactants
  reactant <- unique(unlist(sapply(reactionList,reactants)))
  # Extract all products
  product <- unique(unlist(sapply(reactionList,products)))
  # Possible candidates to be introduced into the system by exchange reactions or by adding more internal reactions.
  orphan <- product[!(product%in%reactant)]
  if (byCompartment == TRUE){
    # Return orphans by compartment
    sapply(compartments(orphan), function(comp){orphan[grep(paste0("\\[",comp,"\\]"),orphan)]}, simplify = FALSE)
  } else {
    # Return all reactants never produced in any reaction.
    return(orphan)
  }
}



