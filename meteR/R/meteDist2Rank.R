
#' @title meteDist2Rank
#'  
#' @description \code{meteESF} calculate the rank distribution of a meteDist object 
#'
#' @details
#' Extracts the predicted rank distribution from a \code{meteDist} object. 
#' This is effectively the quantile function of the distribution. Used, e.g., 
#' in \code{plot.meteDist}
#' 
#' @param x \code{meteDist} object
# @keywords 
#' @export
#' 
#' @examples
#' data(arth)
#' esf1 <- meteESF(spp=arth$spp,
#'                 abund=arth$count,
#'                 power=arth$mass^(.75),
#'                 minE=min(arth$mass^(.75)))
#' sad1 <- sad(esf1)
#' meteDist2Rank(sad1) 
#'                
#' @return A vector of predicted quantiles, typically used to compare against data as in \code{plot.meteDist}
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#  @note other junk to mention
# @seealso add pi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
#  @aliases - a list of additional topic names that will be mapped to this documentation when the user looks them up from the command line.
#  @family - a family name. All functions that have the same family tag will be linked in the documentation.


meteDist2Rank <- function(x) {
  n <- switch(x$type, 
              'sad' = x$state.var['S0'],
              'ssad' = floor(x$state.var['A0']/x$state.var['A']),
              'ipd' = x$state.var['N0'],
              'sipd' = x$state.var['n'])
  x$q(seq(1, 1/n, length=n) - 1/(2*n))
}

