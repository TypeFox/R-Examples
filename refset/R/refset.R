
# TODO/IDEAS:

"
Should I make rs %r% dfr or refset(rs, dfr) work without extra commas?
If so, it has to take the real semantics of dfr, not dfr[] ...

BUG:
weird failures in tests from testthat...
"

#' Create a reference to a subset of an object
#' 
#' Create a refset - a reference to a subset of an object. 
#' When the object changes, the
#' contents of the refset change, and when the refset is changed, the object
#' is changed too.
#' 
#' @param x name of the refset to create, as a bare name or character string
#' @param data the object to refer to
#' @param ... indices to subset with
#' @param drop passed to \code{\link{Extract}}
#' @param dyn.idx update indices dynamically
#' @param read.only create a read-only refset which throws an error if assigned
#'        to
#' @param eval.env environment in which \code{data} and indices will be evaluated
#' @param assign.env environment in which the variable named by \code{x} will be created 
#' 
#' @details
#' There are two ways to call \code{refset}. The two-argument form, e.g.
#' \code{refset(myref, mydata[rows,"column"])}, creates a reference to the 
#' subset of \code{mydata} passed in the second argument. 
#' 
#'The three-or-more argument form acts like the \code{\link{subset}} function:
#'the indices in \code{...} are applied to \code{data}. If \code{data} is a
#'data.frame, then the indices are interpreted within it, so you can refer to
#'columns directly: \code{refset(myref, mydata, a>1 & b<a,)}. Bare column names
#'must be quoted, however.
#'
#' Empty arguments in \code{...} are allowed and are treated as indexing 
#' the whole dimension, just as in \code{\link{Extract}}.
#' 
#' By default, the indices in subset are updated dynamically. 
#' For example, if you call \code{refset(myref, mydata, x >= 3,)} and then
#' set \code{mydata$x <- 3}, the number of rows in \code{myref} will probably
#' increase. To turn this behaviour off and make a reference to a "fixed" 
#' subset of your object, use \code{dyn.idx=FALSE}.
#' 
#' \code{\%r\%} is an infix version of the two-argument form.
#' 
#' @return \code{refset} returns \code{NULL}, but the \code{x} argument 
#' will be assigned to
#' in the calling environment (or in \code{env}, if it is specified). 
#' \code{x} will have an attribute \code{".refset."}.
#' 
#' @seealso
#' Refsets are implemented using \code{makeActiveBinding}.
#' 
#' @aliases `%r%`
#' @family parcel
#' 
#' @examples
#' dfr <- data.frame(a=1:4, b=1:4)
#' ss <- dfr[1:2,]
#' refset(rs, dfr[1:2,])
#' dfr$a <- 4:1
#' ss # 1:2
#' rs # 4:3
#' 
#' # same:
#' refset(rs, dfr, 1:2, )
#' 
#' # same:
#' rs %r% dfr[1:2,]
#' 
#' vec <- 1:10
#' refset(middle, vec[4:6])
#' vec[4:6] <- NA
#' middle 
#' middle <- 4:6 + 100
#' vec
#' 
#' # dynamic versus static indices:
#' dfr <- data.frame(a=rnorm(100), b=rnorm(100))
#' refset(ss, dfr, a>1,)
#' refset(ss.static, dfr, a>1,, dyn.idx=FALSE)
#' nrow(ss) == nrow(ss.static)
#' dfr$a <- dfr$a + 2 * dfr$b
#' 
#' 
#' precious.data <- rnorm(100)
#' refset(big, precious.data, precious.data>1, read.only=TRUE)
#' big
#' \dontrun{
#' big <- big * 2 # throws an error
#' }
#'
#' # Using refset with other functions:
#' # dynamically updated calculated column
#' dfr <- data.frame(a=rnorm(10), b=rnorm(10))
#' refset(rs, transform(dfr, x=a+2*b+rnorm(10)))
#' rs
#' rs # different
#' 
#' # Non-readonly refset with other functions. Works but gives a warning:
#' \dontrun{
#' vec <- 1:5
#' refset(ssv, names(vec), read.only=FALSE)
#' ssv <- LETTERS[1:5]
#' vec
#' }
#' 
## the below is not yet true :-( 
## Works nicely with magrittr or dplyr:
# \dontrun{
# a %>% refset(mydata[1:3,1:4])
# }
#
#' @export
refset <- function(x, data, ..., drop=TRUE, dyn.idx=TRUE, read.only=FALSE,
      eval.env=parent.frame(), assign.env=parent.frame()) {
  env <- eval.env
  if (missing(...)) {
    ssarg <- as.character(substitute(data)[[1]])
    dots <- as.list(substitute(data)[-1:-2])
    data <- substitute(data)[[2]]
    if ("drop" %in% names(dots)) {
      drop <- eval(dots$drop, env) # should it be in env? think so
      dots$drop <- NULL
    }
  } else {
    ssarg <- "["
    dots <- match.call(expand.dots=FALSE)$...
    data <- substitute(data)
  }
  rdata <- eval(data, env)
  if (is.name(substitute(x))) x <- deparse(substitute(x))
  assignargs <- list("["="[<-", "[["="[[<-", "$"="$<-")
  assignarg <- assignargs[[ssarg]]
  if (is.null(assignarg)) {
    if (! missing(read.only) && ! read.only) {
      assignarg <- paste0(ssarg, "<-") 
      warning("Using read.only=FALSE with a non-standard argument of ", 
            sQuote(assignarg), ", assigning to ", sQuote(x), 
            " may cause unexpected behaviour.")
    } else {
      read.only <- TRUE
    }
  }
  missings <- sapply(dots, function(x) is.symbol(x) && identical(
        as.character(x), ""))
  if (length(missings)==0) missings <- logical(0)
  # use argument recycling (of TRUE) as a substitute for non-standard eval:
  dots[missings] <- if (dyn.idx) TRUE else lapply(dim(rdata)[missings], seq, 
        from=1)
  isdfr <- is.data.frame(rdata)
  idxval <- if (dyn.idx) {
    eval(substitute(substitute(dots)), env) 
  } else if (isdfr) {
    lapply(dots, eval, envir=eval(data, env), enclos=env)
  } else {
    lapply(dots, eval, envir=env)    
  }
  
  f <- function(v) {
    if (dyn.idx && isdfr && ssarg != "$") idxval <- lapply(idxval, eval, 
          eval(data, env), env) 
    if (missing(v)) {
      args <- c(data, idxval)
      if (length(idxval) > 1) args <- c(args, drop=drop)
      do.call(ssarg, args, envir=env)
    } else {
      if (read.only) stop("Tried to assign to a readonly refset")
      do.call("<-", list(data, do.call(assignarg, c(data, idxval, 
            list(value=v)), envir=env)), envir=env)
    }
  }
  
  if (exists(x, where=assign.env, inherits=FALSE)) rm(list=x, pos=env)
  makeActiveBinding(x, f, assign.env)
}

#' @export
#' @rdname refset
# `%r%` <- refset
# the broken way to do it, helps pass cran checks, needs new argument to refset
`%r%` <- function(x, data) {
  mc <- match.call()
  mc[[1]] <- quote(refset)
  eval(mc, parent.frame())
}

#' Wrap an expression and its environment into a parcel.
#' 
#' Refsets (and other active bindings) cannot be passed as function
#' arguments, since doing so makes a copy. \code{wrap} allows you to pass
#' arbitrary expressions between functions and records where they are
#' ultimately evaluated. 
#' 
#' @param expr an R expression
#' @param env environment in which \code{expr} is to be evaluated
#' @return
#' An object of class 'parcel', with components \code{expr} and \code{env}.
#' 
#' @family wrapping functions
#' 
#' @examples
#' dfr <- data.frame(a=1:4, b=1:4)
#' rs %r% dfr[1:2,]
#' parcel <- wrap(rs)
#' f <- function (parcel) contents(parcel) <- contents(parcel)*2
#' f(parcel)
#' contents(parcel)
#' dfr
#' 
#' parcel <- wrap(x^2) # non-refset use
#' x <- 3           
#' f <- function(parcel) {x <- 10; contents(parcel)}
#' f(parcel)
#' 
#' @export
wrap <- function(expr, env=parent.frame()) {
  stopifnot(is.environment(env))
  parcel <- new.env(parent=env) # is parent=env necessary?
  parcel$env <- env
  expr <- match.call()$expr
  parcel$expr <- expr
  class(parcel) <- c("parcel", class(parcel))
  parcel
}

#' Convenience function to create a parcel containing a refset.
#' 
#' \code{wrapset} calls \code{\link{refset}} on its arguments and
#' returns the resulting active binding in a \code{parcel} object
#' for passing around.
#' 
#' @param data,... passed to \code{\link{refset}} 
#' @param env passed to \code{\link{refset}} as argument \code{eval.env}
#' 
#' @return A \code{parcel} object.
#' @family wrapping functions
#' 
#' @examples
#' dfr <- data.frame(a=1:5, b=1:5)
#' parcel <- wrapset(dfr, a<3, , drop=FALSE)
#' contents(parcel)
#' 
#' @export
wrapset <- function(data, ..., env=parent.frame()) {
  stopifnot(is.environment(env))
  parcel <- new.env(parent=env) # is parent=env necessary?
  parcel$env <- env
  mc <- match.call(expand.dots=TRUE)
  mc$env <- NULL
  mc$x <- quote(expr)
  mc$eval.env=env
  mc$assign.env=parcel
  mc[[1]] <- quote(refset)
  eval(mc)
  class(parcel) <- c("parcel", class(parcel))
  parcel
}


#' Checks whether an object is a parcel
#' @export
#' @param x an object to examine
#' @return \code{TRUE} or \code{FALSE}.
#' @family wrapping functions
is.parcel <- function(x) inherits(x, "parcel")


#' Returns or changes parcel contents
#' 
#' \code{contents} returns the value of the parcel contents by evaluating
#' the expression in the parcel. \code{contents<-} attempts to assign
#' to the expression, which will only work if the expression is appropriate, e.g.
#' a refset.
#' 
#' @export
#' @param parcel an object of class 'parcel'
#' @param value a value to assign
#' @return The result of evaluating the expression stored in the parcel. 
#' For \code{contents<-}, the parcel itself.
#' 
#' \code{contents<-} will only work if the expression wrapped in the
#' parcel can accept assignments.
#' 
#' @examples
#' pcl <- wrap(x^2)
#' x <- 2
#' contents(pcl)
#' x <- 3
#' contents(pcl)
#' \dontrun{
#' contents(pcl) <- 4 # fails
#' }
#' p2 <- wrap(names(x))
#' contents(p2) <- "named"
#' x
#' 
#' @family wrapping functions
contents <- function(parcel) {
  stopifnot(is.parcel(parcel))
  eval(parcel$expr, parcel$env)
}

#' @export
#' @rdname contents
`contents<-` <- function(parcel, value) {
  stopifnot(is.parcel(parcel))
  expr2 <- bquote(expr <- .(value))
  eval(expr2, parcel)
  parcel
}

#' Unwrap contents of a parcel into a new variable
#' 
#' \code{unwrap_as} creates a new variable which, when evaluated,
#' calls \code{\link{contents}} to return the parcel contents.
#' 
#' @param x name of the variable to bind to
#' @param parcel an object of class 'parcel'
#' @param env environment to assign the variable into
#' @export
#' @family wrapping functions
#' @examples
#' vec <- 1:10
#' parcel <- wrapset(vec, vec > 3) 
#' unwrap_as(y, parcel)
#' y
unwrap_as <- function(x, parcel, env=parent.frame()) {
  stopifnot(is.parcel(parcel))
  x <- deparse(substitute(x))
  f <- function(val) {
    if (missing(val)) {
      contents(parcel) 
    } else {
      contents(parcel) <- val  
    }
  }
  if (exists(x, where=env, inherits=FALSE)) rm(list=x, pos=env)
  makeActiveBinding(x, f, env=env) 
}
