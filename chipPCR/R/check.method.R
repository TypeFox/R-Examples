#check if chosen method is present in method.names vector, use smart autocompletition

check.method <- function(method.names, input.method) {
  res <- method.names[grepl(tolower(input.method), method.names)]
  
  if (length(res) == 0)
    stop("Name ", input.method, " cannot be associated with any available method.")
  
  if (length(res) > 1)
    stop("Name ", input.method, " is too ambiguous. Rerun with more precise name.")
  
  res
}