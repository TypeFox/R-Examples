effect.matrix.factor <-
function (factors, assign = array(dim = 0), index_factor = NA, numeric_index = array(dim = 0)){
  if (length(assign) != 0){
    index <- which(assign == index_factor)
    level <- length(index) + 1
    effect_matrix <- matrix(0, nrow = level, ncol = length(assign))
    effect_matrix[,index] <- contr.sum(level)
    effect_matrix[,1] <- 1 # grand mean included
    #if (length(numeric_index) > 0) effect_matrix[, numeric_index] <- 1 # consider the covariates effect
    attr(effect_matrix, 'levels') <- factors
  }
  return(effect_matrix)
}
