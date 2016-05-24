#' @title Print summaries of \code{meteDist} objects
#'
#' @description
#' S3 method for class \code{meteDist}
#'
#' @details
#' Prints state variables and lagrange multipliers
#' 
#' @param x a \code{meteDist} object (e.g. from \code{ipd.mete} or \code{sad.mete})
#' @param ... arguments to be passed
# @keywords manip
#' @export
#' 
#' @examples
#' data(arth)
#' esf1 <- meteESF(spp=arth$spp,
#'                 abund=arth$count,
#'                 power=arth$mass^(.75),
#'                 minE=min(arth$mass^(.75)))
#' ipd1 <- ipd(esf1)
#' ipd1
#' @return The \code{meteDist} object is returned invisibly
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# @seealso sad.mete, metePsi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

## TODO? add reporting of a data summary (e.g. from str(x) and the parameters of the model (betat and gamma))

print.meteDist <- function(x,...) {
  cat(switch(x$type,
             'sad' = 'Species abundance distribution',
             'ipd' = 'Individual metabolic rate distribution',
             'spd' = 'Species metabolic rate distribuiton',
             'sipd' = 'Species level metabolic rate distribuiton',
             'ssad' = 'Spatial species abundance distribution'), 
      sprintf('predicted using %s', ifelse(is.null(x$data), 
                                           'state variables only', 
                                           'raw data')),
              '\nwith parameters: \n')
  
  print(round(x$state.var, 3))
  print(round(x$La, 5))
  invisible(x)
}
