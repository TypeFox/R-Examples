#' Parse Trellis Formula for \code{strucplot}
#'
#' This is a wrapper for \code{latticeParseFormula} that allows the option of using
#' a "." for the conditioning variables (i.e. after "|") instead of explicitly
#' writing them out. When so used, it means "all variables in the data argument
#' \emph{except} those already used to the left of the |". See the Help for
#' \code{\link{strucplot}} for examples.
#'
#' Note that this is a convenience option only; the
#' conditioning can always be explicitly given. Also note that the two options
#' cannot be combined: either a "." and \emph{only} a "." must be used or
#' \emph{all} the conditioning variables must be written out.
#'
#' @seealso \code{\link[lattice]{latticeParseFormula}}
#'  \code{\link[lattice]{xyplot}}
#'
#' @param form The \code{strucplot} formula to be parsed.
#' @param data An optional data frame containing values for any variables in
#'  the formula. Default = list(), which means that all variables will be looked
#'  up in the formula's environment.
#'
#' @return Same as \code{latticeParseFormula} with an additional "form"
#'  attribute that is the formula used with all conditioning explicitly
#'  given.
#'
#' @examples
#'  exdat <- data.frame(x = 1:5, alongname = sample( letters[1:3],5, rep=TRUE),
#'    butalongername = sample(LETTERS[1:2],5, rep = TRUE))
#'    y <- runif(5)
#'
#'  strucParseFormula (y ~ x| alongname*butalongername, data = exdat)
#'
#'  # The same
#'  strucParseFormula (y ~ x|., data = exdat)
#'
#'  # The 'data' argument is required with '.'
#'  \dontrun{
#'  strucParseFormula (y ~ x|.)
#'  }
#'
strucParseFormula <- function(form, data = list())
  ## parses lattice formulas that can have a '.'
  ## for the conditioning variables in the fomula
{
  badform <- "Bad formula argument"
  if(!inherits(form,"formula"))stop(badform)
  an <- tryCatch(all.names(form),error = function(e)stop(badform))
  if(!("|" %in% an)) stop ("No conditioning variables")
  av <- all.vars(form, unique = FALSE)
  if(any(duplicated(av)))stop("Duplicated variable names not allowed")
  if(any(!(av %in% c(names(data),names(environment(form)),"."))))
    stop("Cannot find variables in formula")
  lenform <- length(form)
  if(!lenform %in% c(2L,3L)) stop(badform)
  else {
    len <- length(an)
    if (an[len] == "." ){
      if(!length(data)) stop("'data' argument required for '.' conditioning variables")
      if(!is.data.frame(data))
         data <- tryCatch(data.frame(data), error=function(e)(
                 stop("Unable to coerce 'data' to a data frame")))
      rightvars <-setdiff(names(data), av[-length(av)])
      if(!length(rightvars)) stop("No conditioning variables")
      else
        form <- eval(do.call(substitute,list(form, env =list(
                      . = parse(text = paste0(rightvars,collapse="*"))[[1]]))))
    }
  }
  out <- latticeParseFormula(form,data)
  cond <- lapply(out$condition,function(x)if(is.shingle(x))ordered(x) else x)
  isfac <- sapply(cond,is.factor)
  if(!all(isfac)){
    bad <- names(cond)[!isfac]
    stop(sprintf("Could not coerce any of %s to factors",
                 paste(bad,collapse = ", ")))
  }
  else {
    out$condition <- cond
    attr(out,"form")<- form
    out
  }
}


