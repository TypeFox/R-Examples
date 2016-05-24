# toChEBI
# Daniel Camilo Osorio
# Bioinformatics and Systems Biology Lab      | Universidad Nacional de Colombia
# Experimental and Computational Biochemistry | Pontificia Universidad Javeriana

# Converts metabolite names to ChEBI ids in a stoichiometric reaction
toChEBI <- function(reaction, formula = FALSE) {
  # Evaluates reversibility
  reversible <- grepl("<=>",reaction)
  if (reversible) {
    reaction<-unlist(strsplit(reaction,"<=>",fixed = TRUE))
    reaction <- gsub("^[[:blank:]]*","",reaction)
    reaction <- gsub("[[:blank:]]*$","",reaction)
  } else {
    reaction<-unlist(strsplit(reaction,"=>",fixed = TRUE))
    reaction <- gsub("^[[:blank:]]*","",reaction)
    reaction <- gsub("[[:blank:]]*$","",reaction)
  }
  reactant <- reaction[1]
  product <- reaction[2]
  reactant <- unlist(strsplit(reactant,"[[:blank:]]\\+[[:blank:]]"))
  product <- unlist(strsplit(product,"[[:blank:]]\\+[[:blank:]]"))
  r_coef <- suppressWarnings(as.numeric(sapply(reactant, .coeficients)))
  r_coef[is.na(r_coef)] <- 1
  p_coef <- suppressWarnings(as.numeric(sapply(product, .coeficients)))
  p_coef[is.na(p_coef)] <- 1
  
  if (formula == FALSE){
    reactant <- mapply(function(coef,met){paste(coef,chebi.id(.metname(met,rm.coef = TRUE)),collapse =" ")}, coef=r_coef,met=reactant)
    product <- mapply(function(coef,met){paste(coef,chebi.id(.metname(met,rm.coef = TRUE)),collapse =" ")}, coef=p_coef,met=product)
  } else {
    reactant <- mapply(function(coef,met){paste(coef,chebi.formula(.metname(met,rm.coef = TRUE)),collapse =" ")}, coef=r_coef,met=reactant)
    product <- mapply(function(coef,met){paste(coef,chebi.formula(.metname(met,rm.coef = TRUE)),collapse =" ")}, coef=p_coef,met=product)
  }
  
  if (reversible){
    paste(paste0(reactant,collapse = " + "),paste0(product, collapse = " + "), sep = " <=> ")
  } else {
    paste(paste0(reactant,collapse = " + "),paste0(product, collapse = " + "), sep = " => ")
  }
}

