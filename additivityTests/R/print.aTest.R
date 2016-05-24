#' Print aTest object
#' 
#' Print the results of additivity test.
#'
#' @param x aTest object
#' @param how many significant digits are to be used
#' 
#' @keywords internal
#'
#' @export
#' 
#' @method print aTest
#' 
#' @return \code{NULL}
#' 
#' @examples
#' data(Boik)
#' t <- tukey.test(Boik)
#' print(t) 

print.aTest<-function(x, digits = max(3, getOption("digits") - 3),...)
{
cat(paste('\n',x$name,' on ',format(x$alpha*100, digits = digits),'% alpha-level:\n\n',sep=""));
cat(paste('Test statistic:',format(x$stat, digits = digits),"\n"));
cat(paste('Critival value:',format(x$critical.value, digits = digits),"\n"));
if (x$result) cat("The additivity hypothesis was rejected.\n\n")
  else cat("The additivity hypothesis cannot be rejected.\n\n")
}
