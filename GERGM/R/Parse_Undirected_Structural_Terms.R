parse_undirected_structural_terms <- function(
  formula,
  possible_structural_terms,
  possible_structural_terms_undirected){

  formula<- as.formula(formula)
  parsed <- deparse(formula)
  if(length(parsed) > 1){
    parsed <- paste0(parsed, collapse = " ")
  }

  for(i in 1:length(possible_structural_terms)){
    if(grepl(possible_structural_terms[i], parsed)){
      valid_term = FALSE
      for(j in 1:length(possible_structural_terms_undirected)){
        if(grepl(possible_structural_terms[i], paste0(possible_structural_terms_undirected,collapse = " "))){
          valid_term = TRUE
        }
      }

      #deal with the case where we have to transform
      if(!valid_term){
        stop(paste("If using an undirected network, you may only specify the following structural terms: edges, ttriads, twostar. You specified the term",possible_structural_terms[i], "."))
      }
    }
  }

  # hand rolled replacements
  parsed <- stringr::str_replace_all(parsed, "twostar","in2star")
  cat("Internally altered formula:", parsed, "\n")
  formula <- as.formula(parsed)
  return(formula)
}
