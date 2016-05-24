check.input <-
function(original.input,col.input=NULL,usage="'input'") {
  col.location <- -9
  if (is.data.frame(original.input)) {
    if(!check.value(col.input,"str_int")) {
      cat(sprintf("Error in %s: Column must be an integer or a character string\n",usage))
      return()
    }
    if (!is.null(col.input) && is.character(col.input)) {
      col.location <- which(names(original.input)==col.input)
      if (length(col.location)!=1) {
        cat(sprintf("Error in %s: Column not found\n",usage))
        return()
      }
    } else if (!is.null(col.input) && is.numeric(col.input)) {
      if (col.input < 1 || col.input > dim(original.input)[2]) {
        cat(sprintf("Error in %s: Column index out of bounds\n",usage))
        return()
      } else {
        col.location <- col.input
      }
    } else if (is.null(col.input) && dim(original.input)[2]==1) {
      col.location <- 1
    } else {
      cat(sprintf("\nPlease specify the index or name of the column in %s containing names\n",usage))
      return()
    }
    vec.names <- original.input[,col.location]
  } else if (check.value(original.input,"vector")) {
    col.location <- 0
    vec.names <- original.input
  } else {
    cat(sprintf("\nPlease enter either a data frame or a vector for %s\n",usage))
    return()
  }  
  if (!is.character(vec.names) && !is.factor(vec.names)) {
    cat(sprintf("Error: Specified vector in %s not alphabetic\n",usage))    
    return()
  }
  return(col.location)  
}
