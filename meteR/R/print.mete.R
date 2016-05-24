#' @title print.meteESF
#'  
#' @description \code{print.meteESF} prints an object of class \code{meteESF}
#'
#' @details
#' See Examples
#' 
#' @param x an object of class \code{meteESF}
#' @param ... arguments to be passed
#' @keywords lagrange multiplier, METE, MaxEnt, ecosystem structure function
#' @export
#' 
#' @examples
#' data(arth)
#' esf1 <- meteESF(spp=arth$spp,
#'                 abund=arth$count,
#'                 power=arth$mass^(.75),
#'                 minE=min(arth$mass^(.75)))
#' print(esf1)
#' esf1 # alternatively...
#' 
#' @return \code{x} silently
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#  @note other junk to mention
# @seealso 
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
#  @aliases - a list of additional topic names that will be mapped to this documentation when the user looks them up from the command line.
#  @family - a family name. All functions that have the same family tag will be linked in the documentation.

##  print and summary methods for METE object

print.meteESF <- function(x,...) {
  cat("METE object with state variables:\n")
  print(x$state.var)
  
  cat("\n");cat("with Lagrange multipliers:\n")
  print(x$La)
  
  #   cat("\n");cat("predicting the distributions:\n")
  #   cat(names(x)[-1]);cat("\n")
}

# summary.METE <- function(x) {

# }

#===========================================================================
# @title print.metePi
#  
# @description \code{print.metePi} prints an object of class \code{metePi}
# 
# @details
# 
# @param x an object of class \code{metePi}
# @param ... arguments to be passed
# @keywords lagrange multiplier, METE, MaxEnt, ecosystem structure function
# @export
# 
# @examples
# 
# @return \code{x} silently
# 
# @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#  @note other junk to mention
# @seealso 
# @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
#  @aliases - a list of additional topic names that will be mapped to this documentation when the user looks them up from the command line.
#  @family - a family name. All functions that have the same family tag will be linked in the documentation.

##  print and summary methods for METE object

## now depricated because output from `metePi' inherits from class `meteESF'
# print.metePi <- function(x,...) {
  # cat("METE object with state variables:\n")
  # print(x$state.var)
  
  # cat("\n");cat("with Lagrange multipliers:\n")
  # print(x$La)
# }
