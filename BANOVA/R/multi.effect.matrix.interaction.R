multi.effect.matrix.interaction <-
function (n_choice, interaction_factors = list(), assign = array(dim = 0), 
                                       main_eff_index = array(dim = 0), index_inter_factor = NA, 
                                       numeric_index = array(dim = 0)){
  result <- list()
  # two way interactions
  if (length(interaction_factors) == 2){
    factor1 <- interaction_factors[[1]]
    factor2 <- interaction_factors[[2]]
    n1 <- length(factor1)
    n2 <- length(factor2)
    temp1 <- array(NA,dim = c(n1*n2,1))
    temp2 <- array(NA,dim = c(n1*n2,1))
    for (i in 1:n1)
      for (j in 1: n2){
        temp1[(i-1)*n2 + j] <- factor1[i]
        temp2[(i-1)*n2 + j] <- factor2[j]
      }
    tempdata <- data.frame(cbind(temp1,temp2)) # everything converted to 1,2,3...
    tempdata[,1] <- as.factor(tempdata[,1])
    tempdata[,2] <- as.factor(tempdata[,2])
    tempmatrix <- model.matrix(~ temp1 * temp2, data = tempdata)
    tempIntermatrix <- tempmatrix[,-(1:(n1-1+n2-1+1))] # intercept:first column is excluded
    tempMainmatrix <- tempmatrix[,2:(n1-1+n2-1+1)]
    index <- which(assign == index_inter_factor)
    for (n_c in 1: n_choice){
      effect_matrix <- matrix(0, nrow = n1*n2, ncol = length(assign))
      effect_matrix[,index] <- tempIntermatrix
      if (n_c != 1)
        effect_matrix[,n_c - 1] <- 1 # grand mean included
      #if (length(numeric_index) > 0) effect_matrix[, numeric_index] <- 1 # consider the covariates effect
      #include main effects
      effect_matrix[, c(which(assign == main_eff_index[1]), which(assign == main_eff_index[2]))] <- tempMainmatrix
      attr(effect_matrix, 'levels') <- cbind(temp1,temp2)
      result[[n_c]] <- effect_matrix
    }
  }else{
    result <- NA
  }
  return(result)
}
