#' \code{rodeo} Reference Class
#'
#' This documents the \code{rodeo} reference class to represent an ODE-based
#' model. See the \code{\link{rodeo-package}} main page or type
#' \code{help(package="rodeo")} for an introduction to the package of
#' the same name.
#'
#' @name rodeo-class
#' @name aliases rodeo-class
#'
#' @field .pros A data frame with fields 'name', 'unit', 'description', and
#'   'expression' defining the process rates.
#' @field .stoi A data frame with fields 'variable', 'process', and 'expression'
#'   reprenting the stoichiometry matrix in data base format.
#' @field .vars A data frame with fields 'name', 'unit', 'description' declaring
#'   the state variables of the model. The declared names become valid
#'   identifiers to be used in the expression fields of \code{.pros} or \code{.stoi}.
#' @field .pars A data frame of the same structure as \code{vars} declaring the
#'   parameters of the model. The declared names become valid
#'   identifiers to be used in the expression fields of \code{.pros} or \code{.stoi}.
#' @field .funs A data frame of the same structure as \code{vars} declaring any
#'   functions referenced in the expression fields of \code{.pros} or \code{.stoi}.
#'
#' @seealso See the \code{\link{rodeo-package}} main page or type
#'   \code{help(package="rodeo")} to find the documentation of any non-class
#'   methods contained in the \code{rodeo} package.
#' 
#' @examples
#' data(exampleIdentifiers, exampleProcesses, exampleStoichiometry)
#' model= new("rodeo",
#'   vars=subset(exampleIdentifiers, type=="v"),
#'   pars=subset(exampleIdentifiers, type=="p"),
#'   funs=subset(exampleIdentifiers, type=="f"),
#'   pros=exampleProcesses, stoi=exampleStoichiometry
#' )
#' model$show()
#'
#' @export

rodeo= setRefClass(
  Class = "rodeo",
  fields = c(
    .vars="data.frame",
    .pars="data.frame",
    .funs="data.frame",
    .pros="data.frame",
    .stoi="data.frame")
)

