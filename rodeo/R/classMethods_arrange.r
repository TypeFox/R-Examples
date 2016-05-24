
# Converts a list into vector after sorting/selecting columns and data recycling
arrangeGeneric= function(x, itemNames, sep=".", warnUnused=TRUE) {
  # Check inputs
  if (!is.list(x))
    stop("actual argument must be a list")
  if (length(x) == 0)
    stop("input list has length zero")
  bad= itemNames[!(itemNames %in% names(x))]
  if (length(bad) > 0)
    stop(paste0("input list does not provide data for the following",
      " item(s): '",paste(bad,collapse="', '"),"'"))
  if (warnUnused) {
    bad= names(x)[!(names(x) %in% itemNames)]
    if (length(bad) > 0)
      warning(paste0("the following element(s) from the input list were ignored: '",
        paste(bad,collapse="', '"),"'"))
  }
  # Transform into data frame to get the same number of levels for all variables
  sizes= unique(sapply(x, length))
  if ((length(sizes) > 2) || ((length(sizes) == 2) && (min(sizes) != 1)))
    stop("elements of input list must either be scalars or vectors of a common length")
  x= data.frame(x)
  nr= nrow(x)
  # Sort/drop columns
  x= x[,itemNames]
  # Transform into vector
  x= unlist(x)
  if (nr > 1) {
    names(x)= paste(rep(itemNames, each=nr), rep(1:nr, length(itemNames)), sep=".")
  } else {
    names(x)= itemNames
  }
  # Check type
  if (any(!is.numeric(x)))
    stop("non-numeric data detected in input list")
  # Done
  return(x)
}

#' Assign Values to State Variables
#'
#' Assign values to state variables of a \code{\link{rodeo}}-based model.
#'
#' @name arrangeVars
#'
#' @param x A list of values to be assigned to the model's state variables.
#'   The elements of \code{x} are either numeric scalars (for
#'   single box models) or vectors (for multi-box models). In the latter case,
#'   the length of the vector(s) must be equal to the number of modeled boxes
#'   (e.g. layers in a vertical 1D model).
#'   If scalars and vectors are present, the scalars will be turned into vectors
#'   by replication.
#' @param sep A character string used to create element names for the return
#'   vector in multi-box models (i.e. if some element(s) of \code{x} are
#'   vectors). Element names of the return vector are constructed by appending
#'   the level index (starting from 1) to the name of the respective list
#'   element using the value of \code{sep} as a separator. In the single-box
#'   case, \code{sep} is ignored, and the elements of the return vector carry
#'   the same names as the corresponding list elements.
#' @param warnUnused Logical. If \code{TRUE}, a warning is issued if \code{x}
#'   an element name in \code{x} does not match the name of a state variable.
#'   Setting this to \code{FALSE} allows \code{x} to contain auxiliary data.
#'
#' @return Named numeric vector with elements being ordered for compatibility
#'   with the respective \code{\link{rodeo}} object. The returned vector can be
#'   passed to the \code{y} argument of \code{\link[deSolve]{ode}}, for example.
#'
#' @author \email{david.kneis@@tu-dresden.de}
#'
#' @note No additional notes so far.
#'
#' @seealso Other \code{\link{rodeo}} class methods
#'
#' @examples
#' data(exampleIdentifiers, exampleProcesses, exampleStoichiometry)
#' model= new("rodeo",
#'   vars=subset(exampleIdentifiers, type=="v"),
#'   pars=subset(exampleIdentifiers, type=="p"),
#'   funs=subset(exampleIdentifiers, type=="f"),
#'   pros=exampleProcesses, stoi=exampleStoichiometry
#' )
#' v= model$arrangeVars(list(c_z=1, c_do=9.022, v=1.e6))
#' print(v)

rodeo$methods( arrangeVars= function(x, sep=".", warnUnused=TRUE) {
  "Assign values to state variables. See \\code{\\link{arrangeVars}} for details."
  arrangeGeneric(x=x, itemNames=.self$.vars$name, sep=sep, warnUnused=warnUnused)
})

#' Assign Values to Parameters
#'
#' Assign values to parameters of a \code{\link{rodeo}}-based model.
#'
#' @name arrangePars
#'
#' @param x A list of values to be assigned to the model's parameters.
#'   The elements of \code{x} are either numeric scalars (for
#'   single box models) or vectors (for multi-box models). In the latter case,
#'   the length of the vector(s) must be equal to the number of modeled boxes
#'   (e.g. layers in a vertical 1D model).
#'   If scalars and vectors are present, the scalars will be turned into vectors
#'   by replication.
#' @param sep A character string used to create element names for the return
#'   vector in multi-box models (i.e. if some element(s) of \code{x} are
#'   vectors). Element names of the return vector are constructed by appending
#'   the level index (starting from 1) to the name of the respective list
#'   element using the value of \code{sep} as a separator. In the single-box
#'   case, \code{sep} is ignored, and the elements of the return vector carry
#'   the same names as the corresponding list elements.
#' @param warnUnused Logical. If \code{TRUE}, a warning is issued if \code{x}
#'   an element name in \code{x} does not match the name of a parameter.
#'   Setting this to \code{FALSE} allows \code{x} to contain auxiliary data.
#'
#' @return Named numeric vector with elements being ordered for compatibility
#'   with the respective \code{\link{rodeo}} object. The returned vector can be
#'   passed to the \code{parms} argument of \code{\link[deSolve]{ode}}, for
#'   example.
#'
#' @author \email{david.kneis@@tu-dresden.de}
#'
#' @note No additional notes so far.
#'
#' @seealso Other methods of the \code{\link{rodeo-class}}
#'
#' @examples
#' data(exampleIdentifiers, exampleProcesses, exampleStoichiometry)
#' model= new("rodeo",
#'   vars=subset(exampleIdentifiers, type=="v"),
#'   pars=subset(exampleIdentifiers, type=="p"),
#'   funs=subset(exampleIdentifiers, type=="f"),
#'   pros=exampleProcesses, stoi=exampleStoichiometry
#' )
#' p= model$arrangePars(list(kd=5.78e-7, h_do=0.5, s_do_z=2.76, wind=1, depth=2,
#'   temp=20, q_in=1, q_ex=1))
#' print(p)


rodeo$methods( arrangePars= function(x, sep=".", warnUnused=TRUE) {
  "Assign values to parameters. See \\code{\\link{arrangePars}} for details."
  arrangeGeneric(x=x, itemNames=.self$.pars$name, sep=sep, warnUnused=warnUnused)
})

