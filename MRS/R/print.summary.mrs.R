#' Print summary of a mrs object
#'
#' This function print the summary the output of the \code{\link{mrs}} function.
#' It provides the marginal prior and posterior of the null and the top regions of the representative tree.
#'
#' @param x A \code{summary.mrs} object
#' @param ... Additional print parameters. 
#' @references Soriano J. and Ma L. (2014). Multi-resolution two-sample comparison 
#' through the divide-merge Markov tree. \emph{Preprint}. 
#'  \url{http://arxiv.org/abs/1404.3753}
#' @export
#' @S3method summary mrs
#' @examples
#' set.seed(1) 
#' n = 100
#' p = 2
#' X = matrix(c(runif(p*n/2),rbeta(p*n/2, 1, 4)), nrow=n, byrow=TRUE)
#' G = c(rep(1,n/2), rep(2,n/2))
#' x = mrs(X=X, G=G)
#' fit = summary(x, rho = 0.95, abs_eff = 1)
#' print(fit)
print.summary.mrs <-function(x, ...)
{
  cat("------------------------\n")
  cat("Posterior Null: ", x$Posterior_Null, "\n")
  cat("Prior Null    : ", x$Prior_Null, "\n")
  cat("------------------------\n\n")
  if( x$Num_Regions>0)
  {
    cat("Top differential regions \n")    
    cat("------------------------\n")
    for(i in 1:x$Num_Regions)
    {
      cat("PMAP: ", signif(x$Alt_Prob[i], digits=3), " |\t") 
      cat("Effect Size: ", signif(x$Effect_Size[[i]], digits=3), " |\t")
      cat("Region: ", signif(x$Regions[i,], digits=3), " |\t")
      cat("Direction: ", x$Directions[i], "\n")
      cat("------------------------\n")
    }
  }
  
}