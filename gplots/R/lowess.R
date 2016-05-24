# make stats::lowess into a generic base-function
lowess.default <- function (x, y = NULL,
                            f = 2/3,
                            iter = 3L,
                            delta = 0.01 * diff(range(x)),
                            ...)
  {
    m <- match.call()
    m[[1L]] <- quote(stats::lowess)
    retval <- eval(m, envir=parent.frame())
    class(retval) <- "lowess"
    retval$call <- match.call()
    retval
  }

lowess  <- function(x,...)
  UseMethod("lowess")


"lowess.formula" <-  function (formula,
                               data = parent.frame(),
                               ...,
                               subset,
                               f=2/3,
                               iter=3,
                               delta=.01*diff(range(mf[-response]))
                               )
{
  if (missing(formula) || (length(formula) != 3))
    stop("formula missing or incorrect")

  m <- match.call(expand.dots = FALSE)
  eframe <- parent.frame()
  md <- eval(m$data, eframe)
  if (is.matrix(md))
    m$data <- md <- as.data.frame(data)
  dots <- lapply(m$..., eval, md, eframe)
  nmdots <- names(dots)
  m$...  <- m$f <- m$iter <- m$delta <- NULL
  subset.expr <- m$subset
  m$subset <- NULL
  m <- as.list(m)
  m[[1L]] <- stats::model.frame.default
  m <- as.call(c(m, list(na.action = NULL)))
  mf <- eval(m, eframe)
  if (!missing(subset)) {
    s <- eval(subset.expr, data, eframe)
    l <- nrow(mf)
    dosub <- function(x) if (length(x) == l)
      x[s]
    else x
    dots <- lapply(dots, dosub)
    mf <- mf[s, ]
  }

  mf <- na.omit(mf)

  response <- attr(attr(mf, "terms"), "response" )
  retval <- stats::lowess(mf[[-response]], mf[[response]], f=f, iter=iter, delta=delta)
  class(retval) <- "lowess"
  retval$call <- match.call()

  retval
}
