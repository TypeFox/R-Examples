get.interactions <-
function(mf, ncol_intercepts = 0){
  factor_matrix <- attr(attr(mf,'terms'),'factors')
  factor_order <- attr(attr(mf,'terms'),'order')
  col_index_inter <- which(factor_order >= 2)
  col_chosen_inter <- array(dim = 0) # store interactions without numeric variables
  results <- list()
  if (length(col_index_inter) > 0){
    for (i in 1:length(col_index_inter)){
      index_inter <- which(factor_matrix[,col_index_inter[i]] == 1)
      # exclude the case that the interaction include a numeric variable
      if (sum(attr(attr(mf,'terms'),'dataClasses')[index_inter] != 'factor') == 0){ 
        results[[i]] <- index_inter + ncol_intercepts
        col_chosen_inter <- c(col_chosen_inter, col_index_inter[i] + ncol_intercepts)
      }
    }
  }
  sol <- list()
  sol$results <- results
  sol$index <- col_chosen_inter
  return(sol)
}
