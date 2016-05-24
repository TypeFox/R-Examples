BASIX.unique <- function(matrix){

if(length(rownames(matrix))==0){
   rownames(matrix) <- 1:dim(matrix)[1] 
}


if(is.character(matrix)){
ids    <- .Call("my_unique_C2",matrix,PACKAGE="BASIX")
}else{
ids    <- .Call("my_unique_C",matrix,PACKAGE="BASIX")
}

matrix <- matrix[!ids,]

return(matrix)
}

