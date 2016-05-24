# Method to convert a data.frame with matrix to a standard data.Frame

# Check if there is matrix embended

search_matrix_position <- function(dataframe){
  
  #init
  matrix_position <- integer(0)
  
  # For each variable, search if it's a matrix
  for( variable in 1:ncol(dataframe)){
    classe_of_this_variable <- class(dataframe[,variable])[1]
    if(classe_of_this_variable == "matrix"){
      matrix_position <- append(matrix_position, variable)
    }else{
    }
  }
  
  # Return a logical vector
  return(matrix_position)
}

expand_dfmatrix <- function(dataframe, matrix_var = NA){
  # This problem is complex because there should be more than a matrix
  
  # First version forgot to add the name of the original matrix variable like matrixvar.1
  if(is.na(matrix_var)){
    # if the position of the matrix is not given, find it
    matrix_var <- search_matrix_position(dataframe)
  }
  
  if( length(matrix_var) > 0){
    return(cbind(dataframe[,-c(matrix_var)],dataframe[,matrix_var]))
  }else{
    return(dataframe)
  }
}
