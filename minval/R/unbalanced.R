# unbalanced
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

# Returns the unbalanced reactions from a set of stoichiometric reactions
unbalanced <- function(reaction, show.formulas = FALSE) {
  is.unbalanced <- function(mfreaction){
    if (grepl("<=>",mfreaction)){
      rxn <- unlist(strsplit(mfreaction,"<=>"))
      rxn <- gsub("^[[:blank:]]","",rxn)
      rxn <- gsub("[[:blank:]]$","",rxn)
    } else {
      rxn <- unlist(strsplit(mfreaction,"=>"))
      rxn <- gsub("^[[:blank:]]","",rxn)
      rxn <- gsub("[[:blank:]]$","",rxn)
    }
    reactant <- rxn[1]
    product <- rxn[2]
    reactant <- unlist(strsplit(reactant,"[[:blank:]]\\+[[:blank:]]"))
    reactant <- .atoms(reactant)
    product <- unlist(strsplit(product,"[[:blank:]]\\+[[:blank:]]"))
    product <- .atoms(product)
    return(!identical(.formula2matrix(reactant), .formula2matrix(product)))
  }
  mfreaction <- sapply(reaction,function(reaction){toChEBI(reaction,formula = TRUE)},USE.NAMES = FALSE)
  mb <- sapply(mfreaction, is.unbalanced,USE.NAMES = FALSE)
  if (show.formulas==FALSE){
    return(mb)
  } else{
    ub <- (mb==TRUE)
    cbind(reaction[ub],mfreaction[ub])
  }
}
