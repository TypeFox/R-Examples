trimSpace <- function(x, leading = TRUE, trailing = TRUE, space = "[:space:]"){
  if(leading){
    pattern <- paste("^[", space, "]*", sep = "", collapse = "")
    x <- sub(pattern = pattern, replacement =  "", x = x)
  }
  if(trailing){
    pattern <- paste("[", space, "]*$", sep = "", collapse = "")
    x <- sub(pattern = pattern, replacement =  "", x = x)
  }
  return(x)
}
