print.results <- function(x, ...)
{
  x<-x$output 
  n<-6
  if (dim(x)[2] > n) {
    prettydata <- data.frame(
      head(x)[,c(1:n)], rep('...',dim(head(x))[1])
    )
    colnames(prettydata)[(n+1)] <- '...'
    print (prettydata)
    cat("\nOnly the first few rows and columns are printed.\n")
  } else {
    print(head(x))
    cat("\nOnly the first few rows are printed.\n")
  }
}
