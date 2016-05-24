#Turns the input data that is in the format of a matrix or data frame into a Reference to a 
#java double[][] Object.
#'@import rJava
java_object_from_data <- function(data) {
  #To achieve this, the input data matrix (or data frame) is first turned into a vector column by column. This
  #Vector is then passed into a java function that also needs to know the number of columns of the original
  #matrix to reconstruct it as a double[][]. This is much faster than producing a double[][] with .jarray.
  res <- rJava::.jcall("JavaObjectFromDataConverter",returnSig="[[D",method="matrix_from_array",
                       as.vector(as.matrix(data)),ncol(data),evalArray=F)
  return(res)
}
