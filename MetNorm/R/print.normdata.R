print.normdata <- function(x, ...)
{
  printdata(x$newY)
#  printdata(x$UVcomp)
#  printdata(x$k)
#  print(x$lambda)
#  cat("\nGroups:\n")
#  print(x$groups)
}

printdata<-function(input){
  if (dim(input)[2] > 5) {
    prettydata <- data.frame(head(input)[,c(1:6)], rep('...',dim(head(input))[1]))
    colnames(prettydata)[7] <- '...'
    print (prettydata)
    cat("\nOnly the first few rows and columns are printed.\n")
  } else {
    print(head(input))
    cat("\nOnly the first few rows are printed.\n")
  }
}

