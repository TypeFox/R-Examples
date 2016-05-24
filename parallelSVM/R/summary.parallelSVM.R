#' @export
summary.parallelSVM <- function(object, ...){
  # Calculate average number of support Vectors
  nSV <- c()
  for (i in 1:length(object)){
    nSV <- c(nSV,object[[i]]$tot.nSV)
  }
  tot.nSV <- round(mean(nSV))
  
  # Calculate the division of these support Vectors
  SV <- c()
  for (i in 1:length(object)){
    SV <- c(SV,object[[i]]$nSV)
  }
  tot.SV <- round(mean(SV))
  
  # The normal print elements
  print(object)
  
  # Extra output for summary
  cat("\n\n ( ")
  cat(tot.SV)
  cat(" ")
  cat(tot.nSV - tot.SV)
  cat(" )\n\n\nNumber of classes: ")
  cat(object[[1]]$nclasses)
  cat("\n\nLevels:\n")
  cat(object[[1]]$levels)
  cat("\n\n\n ")
}

