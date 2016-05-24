#' Return the Stoichiometry Matrix
#'
#' Return and optionally evaluate the mathematical expression appearing in the
#' stoichiometry matrix.
#'
#' @name stoichiometry
#'
#' @param values A named numeric vector specifying the values of all state
#'   variables and parameters. For non-autonomous models, there must also be an
#'   element named 'time'. If \code{values} is set to \code{NULL}, the 
#'   mathematical expressions are not evaluated but returned as character
#'   strings.
#'
#' @return A matrix of numeric or character type, depending on the contents of
#'  \code{values}.
#'
#' @note If the stoichiometric factors are mathematical expressions involving
#'   function references, these functions must be defined in R (even if the
#'   numerical computations are based on generated Fortran code).
#'
#' @author \email{david.kneis@@tu-dresden.de}
#'
#' @seealso See other methods of the \code{\link{rodeo-class}} or
#'   \code{\link{plotStoichiometry}} for a graphical representation of the
#'   stoichiometric factors only.
#'
#' @examples
#' data(exampleIdentifiers, exampleProcesses, exampleStoichiometry)
#' model= new("rodeo",
#'   vars=subset(exampleIdentifiers, type=="v"),
#'   pars=subset(exampleIdentifiers, type=="p"),
#'   funs=subset(exampleIdentifiers, type=="f"),
#'   pros=exampleProcesses, stoi=exampleStoichiometry
#' )
#' print(model$stoichiometry())
#' c_z_in= function(time) {0.1}
#' c_do_in= function(time) {8.0}
#' print(model$stoichiometry(c(s_do_z=2.76, c_z=1, c_do=9.022, time=0)))

rodeo$methods( stoichiometry = function(values=NULL) {
  "Returns the stoichiometry matrix, either evaluated (numeric) or not (text).
  See \\code{\\link{stoichiometry}} for details."

  # Build the matrix of expressions
  m= matrix("0", ncol=nrow(.self$.vars), nrow=nrow(.self$.pros))
  colnames(m)= .self$.vars$name
  rownames(m)= .self$.pros$name
  for (i in 1:nrow(.self$.stoi)) {
    m[.self$.stoi$process[i], .self$.stoi$variable[i]]= .self$.stoi$expression[i]
  }

  # Return the matrix of expressions if no values are supplied ...
  if (is.null(values)) {
    return(m)
  # ... or return the numeric matrix otherwise
  } else {
    # Check supplied values
    if (is.null(names(values)) || any(names(values) == ""))
      stop("missing element name(s) in vector 'values'")
    if (any(duplicated(names(values))))
      stop("duplicated element name(s) in vector 'values'")
    if (!all(is.numeric(values)))
      stop("non-numeric element(s) in 'values'")
    if (!all(is.finite(values)))
      stop("non-finite element(s) in 'values'")
    # Create environment holding all data -> required for evaluating expressions
    env= new.env()
    f=tempfile()
    write.table(file=f, x=data.frame(names(values),values,stringsAsFactors=FALSE),
      sep="=", col.names=FALSE, row.names=FALSE, quote=FALSE)
    sys.source(file=f,envir=env)
    # Create numeric matrix
    mnum= matrix(0, ncol=ncol(m), nrow=nrow(m))
    colnames(mnum)= colnames(m)
    rownames(mnum)= rownames(m)
    for (ic in 1:ncol(m)) {
      for (ir in 1:nrow(m)) {
        tryCatch({
          mnum[ir,ic]= eval(parse(text=m[ir,ic]), envir=env)  # evaluated in created env
        }, error= function(e) {
          stop(paste0("failed to compute stoichiometry factor for variable '",
            colnames(m)[ic],"' and process '",rownames(m)[ir],"'; details: ",e))
        })
      }
    }
    return(mnum)
  }
})

