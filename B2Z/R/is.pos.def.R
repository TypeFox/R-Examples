####################################
#This function checks whether a    #
#given matrix is positive definite #
####################################

is.pos.def <- function(mat)
 {
 val = try(chol(mat), silent = TRUE)
 if (class(val) == "try-error") 
 return(FALSE)
 else return(TRUE)
 }