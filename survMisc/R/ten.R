#' @name ten
#' @title \bold{t}ime, \bold{e}vent(s) and \bold{n}umber at risk.
#'
#' @include print.R
#' @include asWide.R
#'
#' @param x
#' For the default method, a \code{numeric} vector indicating an
#' \emph{event} (or status).
#'  \cr
#' Each element indicates whether an event occurred (\code{1}) or
#' not (\code{0}) for an observation.
#'  \cr
#' These are assumed to be ordered by discrete times.
#'  \cr
#' This is similar to the \code{event} argument for \code{Surv}
#' objects.
#'  \cr \cr
#' Methods are available for objects of class
#' \code{Surv}, \code{survfit},
#' \code{coxph} and \code{formula}. 
#' @param abbNames \bold{Abb}reviate names?
#'  \cr
#' If \code{abbNames="TRUE"} (the default), 
#' the covariate groups are referred to by number. 
#'  \cr
#' As the names for each covariate group are made by concatenating
#' the predictor names, the full names can become unwieldly.
#'  \cr
#' If \code{abbNames="FALSE"}, the full names are given.
#'  \cr
#' In either case, the \code{longNames} are given
#' as an \code{attribute} of the returned \code{ten} object.
#' @param contrasts.arg Methods for handling factors.
#'  \cr
#' A \code{list}. The \code{names} are the names of
#' columns of the \code{model.frame} containing
#' \code{factor}s.
#'  \cr
#' The \emph{values} are used as replacement
#' values for the \code{stats::contrasts} replacement function.
#' These should be functions (given as character strings)
#' or numeric matrices. 
#'  \cr
#' This can be passed from
#' \code{survfit}, \code{coxph} and \code{formula} objects to:
#'  \cr
#' ?stats::model.matrix
#' @param call Used to pass the \code{call} from a \code{formula}
#'  to the final \code{ten.data.table} method.
#' @param mm Used to pass the \code{model.matrix} from a \code{formula}
#'  to the final \code{ten.data.table} method.
#' @inheritParams sf.ten
#'
#' @return A \code{data.table} with the additional \code{class}
#' \code{ten}. 
#'  \cr
#' By default, the shape returned is 'long' i.e. there is one row for each unique 
#'  timepoint per covariate group. 
#'  \cr
#' The basic form, for a \code{numeric} or \code{Surv} object, has columns: 
#'  \item{t}{\bold{t}ime.}
#'  \item{e}{number of \bold{e}vents.}
#'  \item{n}{\bold{n}umber at risk.}
#' A \code{survfit}, \code{coxph} or \code{formula} object
#' will have additional columns:
#'  \item{cg}{\bold{c}ovariate \bold{g}roup.
#'   This is formed by combining the variables; these
#'   are separated by a comma ','.}
#'  \item{ncg}{\bold{n}umber at risk, by \bold{c}ovariate \bold{g}roup}
#'  
#' \bold{Special terms}.
#'  \cr \cr
#' The following are considered 'special'
#' terms in a survival model:
#'   \item{strata}{For a stratified model, \code{ten} returns a \code{list} with
#'    one element per strata, which is a \code{ten} object.
#'     \cr
#'    This has the class \code{stratTen}. The name of the 
#'    list elements are those of the strata in the model.}
#'   \item{cluster}{These terms are dropped.}
#'   \item{tt}{The variable is unchanged. That is, time-transform
#'    terms are handled as if the the function
#'    \code{tt(x)} was \code{identity(x)}.}
#' \bold{Attribures}.
#'  \cr
#' The returned object will also have the following \code{attributes}:
#'  \item{shape}{The default is \code{"long"} but 
#'   is changed to \code{"wide"} when \code{asWide} is called on the object.}
#'  \item{abbNames}{Abbreviate names?}
#'  \item{longNames}{A \code{data.table} with two columns, showing the abbrevbiated 
#'   and full names.}
#'  \item{ncg}{Number of covariate groups}
#'  \item{call}{The call used to generate the object}
#'  \item{mm}{The \code{model.matrix} used to generate to 
#'   generate the object, if applicable.}
#' Additional attributes will be added by the following functions:
#'  \cr
#' \code{\link{sf}}
#' \code{\link{ci}}
#'
#' @note 
#' The methods for \code{data.frame} (for a model frame)
#' and \code{data.table} are not typically intended for interactive use.
#'  \cr \cr
#' Currently only binary status and right-censoring
#' are supported. 
#'  \cr \cr
#' In stratified models, only one level of stratification is supported
#' (i.e. strata cannot be 'nested' currently).
#'  \cr \cr
#' Partial matching is available for the
#' following arguments, based on the characters in bold:
#' \itemize{
#'  \item \bold{abb}Names
#'  \item \bold{con}trasts.arg
#' }
#' 
#' @seealso \code{\link{asWide}}
#' @seealso \code{\link{print}}
#'
#' @rdname ten
#' @export
#'
ten <- function(x, ...) UseMethod("ten")
### all are methods ultimately passed to
### ten.data.frame (below)
### except ten.numeric() and ten.Surv()
###----------------------------------------
#'
#' @rdname ten
#' @method ten numeric
#' @aliases ten.numeric
#' @export
#' @examples
#' require("survival")
#' ## binary vector
#' ten(c(1, 0, 1, 0, 1))
#'
ten.numeric <- function(x, ...){
  stopifnot(all(x >= 0 && x <=1))
  res1 <- data.table::data.table(
    "t"=(t <- seq_along(x)),
    "n"=rev(t),
    "e"=x)
  data.table::setattr(res1, "class", c("ten", class(res1)))
  setAttr(res1,
          shape="long",
          abbNames=TRUE,
          ncg=0,
          call=match.call())
  return(res1)
}
#'
#' @rdname ten
#' @method ten Surv
#' @aliases ten.Surv
#' @export
#' @examples
#' ## Surv object
#' df0 <- data.frame(t=c(1, 1, 2, 3, 5, 8, 13, 21),
#'                   e=rep(c(0, 1), 4))
#' s1 <- with(df0, Surv(t, e, type="right"))
#' ten(s1)
#' ## some awkward values
#' suppressWarnings(
#'     s1 <- Surv(time=c(Inf, -1, NaN, NA, 10, 12),
#'                event=c(c(NA, 1, 1, NaN, Inf, 0.75))))
#' ten(s1)
#'
ten.Surv <- function(x, ..., 
                    call=NULL){
  stopifnot(inherits(x, "Surv"))
  stopifnot(attributes(x)$type=="right")
  if(is.null(call)) call <- match.call()
  res1 <- data.table::data.table(unclass(x))
  data.table::setkey(res1, "time")
    res1 <- res1[, list("n"=length(status),
                        "e"=sum(status)),
                 by=sort(time, na.last=TRUE)]
  res1[, "n" := c(sum(n), sum(n) - cumsum(n)[ - length(n)])]
  data.table::setnames(res1, c("t", "n", "e"))
  data.table::setattr(res1, "class", c("ten", class(res1)))
  setAttr(res1,
          shape="long",
          ncg=0,
          call=call)
  return(res1)
}
#'
#' @rdname ten
#' @method ten coxph
#' @aliases ten.coxph
#' @export
#' @examples
#' ## coxph object
#' ## K&M. Section 1.2. Table 1.1, page 2.
#' data("hodg", package="KMsurv")
#' hodg <- data.table::data.table(hodg)
#' data.table::setnames(hodg,
#'                      c(names(hodg)[!names(hodg) %in%
#'                                    c("score", "wtime")],
#'                        "Z1", "Z2"))
#' c1 <- coxph(Surv(time=time, event=delta) ~ Z1 + Z2,
#'             data=hodg[gtype==1 && dtype==1, ])
#' ten(c1)
#' data("bmt", package="KMsurv")
#' ten(c1 <- coxph(Surv(t2, d3) ~ z3*z10, data=bmt))
#' ## T&G. Section 3.2, pg 47.
#' ## stratified model
#' data("pbc", package="survival")
#' c1 <- coxph(Surv(time, status==2) ~ log(bili) + age + strata(edema), data=pbc)
#' ten(c1)
#' 
ten.coxph <- function(x, ...,
                      abbNames=TRUE,
                      contrasts.arg=NULL){
  partMatch(env1=environment(), ...)
  x$call$formula <- stats::terms(
      x=stats::formula(x$call),
      specials=c("strata", "cluster", "tt"))
  mode(x$call) <- "list"
  length(x$call) <- 3
  mode(x$call) <- "call"
  call1 <- x$call
  x$call$drop.unused.levels <- TRUE
  x$call[[1]] <- as.name("model.frame")
  ## model.frame
  xMF1 <- eval(x$call, parent.frame())
  ten(x=xMF1,
      abbNames=abbNames,
      contrasts.arg=contrasts.arg,
      call=call1)
}
#' @rdname ten
#' @aliases ten.survfit
#' @method ten survfit
#' @export
#' @examples
#' ## K&M. Example 7.2, pg 210.
#' data("kidney", package="KMsurv")
#' with(kidney[kidney$type==2, ], ten(Surv(time=time, event=delta)))
#' s1 <- survfit(Surv(time=time, event=delta) ~ type, data=kidney)
#' ten(s1)[e > 0, ]
#'
ten.survfit <- function(x, ..., 
                        abbNames=TRUE,
                        contrasts.arg=NULL){
  partMatch(env1=environment(), ...)
  x$call$formula <- stats::terms(
    x=stats::formula(x$call),
    specials=c("strata", "cluster", "tt"))
  mode(x$call) <- "list"
  length(x$call) <- 3
  mode(x$call) <- "call"
  call1 <- x$call
  x$call$drop.unused.levels <- TRUE
  x$call[[1]] <- as.name("model.frame")
  xMF1 <- eval(x$call, parent.frame())
  ten(x=xMF1,
      abbNames=abbNames,
      contrasts.arg=contrasts.arg,
      call=call1)
}
#'
#' @rdname ten
#' @method ten formula
#' @aliases ten.formula
#' @export
#' @examples
#' ## A null model is passed to ten.Surv
#' (t1 <- with(kidney, ten(Surv(time=time, event=delta) ~ 0)))
#' ## but the original call is preserved
#' attr(t1, "call")
#' ## survival::survfit doesn't accept interaction terms...
#' \dontrun{
#'     s1 <- survfit(Surv(t2, d3) ~ z3*z10, data=bmt)}
#' ## but ten.formula does:
#' ten(Surv(time=t2, event=d3) ~ z3*z10, data=bmt)
#' ## the same is true for the '.' (dot operator) in formulas
#' (t1 <- ten(Surv(time=t2, event=d3) ~ ., data=bmt))
#' ## impractical long names stored as an attribute
#' attr(t1, "longNames")
#'
ten.formula <- function(x, ...,
                        abbNames=TRUE,
                        contrasts.arg=NULL){
  partMatch(env1=environment(), ...)
  stopifnot(inherits(x, "formula"))
  ## based on code from stats::lm()
  mc1 <- match.call()
  names(mc1)[names(mc1)=="x"] <- "formula"
  mc1 <- mc1[c(1L, match(c("formula", "data"), names(mc1), 0L))]
  call1 <- mc1
  mc1$drop.unused.levels <- TRUE
  mc1[[1L]] <- as.name("model.frame")
  mf1 <- eval(mc1, parent.frame())
  ten(x=mf1,
      abbNames=abbNames,
      contrasts.arg=contrasts.arg,
      call=call1)
}
#'
#' @rdname ten
#' @method ten data.frame
#' @aliases ten.data.frame
#' @export
#' @examples
#' ## not typically intended to be called directly
#' mf1 <- model.frame(Surv(time, status==2) ~ age + strata(edema) + strata(spiders), pbc, 
#'                    drop.unused.levels = TRUE)
#' ten(mf1)
#' 
ten.data.frame <- function(x, ...,
                           abbNames=TRUE,
                           contrasts.arg=NULL,
                           call=NULL){
    stopifnot(survival::is.Surv(x[[1]]))
    stopifnot(attr(x[[1]], "type") == "right")
    partMatch(env1=environment(), ...)
    if (stats::is.empty.model(stats::terms(x))) {
        ## extract Surv object
        return(ten(x[[1]],
                   call=call))
    }
    ## names of strata
    xNS1 <- grep("^strata\\(.*\\)", names(x))
    ## data.table from x
    xDT <- data.table::as.data.table(
        stats::model.matrix(stats::terms(x),
                            x,
                            contrasts.arg=contrasts.arg))
    xDT[, c("time", "status") :=
        lapply(1:2L, function(i) stats::model.response(x)[, i])]
    ## names of clusters
    xNC1 <- grep("^cluster\\(.*\\)", names(x))
    if (any(xNC1)) {
        ## drop cluster terms
        data.table::set(xDT,
                        j=grep("^cluster\\(.*\\)", names(xDT)),
                        value=NULL)
    }
    if (any(xNS1)) {
        ## strata numbers
        xDTnS1 <- grep("^strata\\(.*\\)", names(xDT))
#        xDTSn1 <- grep("^strata\\(.*\\)", names(xDT), value=TRUE)
        ## separate table only for strata
        xDTstr <- xDT[, .SD, .SDcols=xDTnS1]
        data.table::set(xDT, j=xDTnS1, value=NULL)
        setnames(xDTstr, sub("^strata\\(.+\\)", "", names(xDTstr)))
        c1 <- colnames(xDTstr)
        xDTstr[, (c1) := lapply(.SD, as.logical), .SDcols=seq.int(ncol(xDTstr))]
        collapseDT(xDTstr, except=NA, nName="strat")
        collapseDT(xDT,
                   except=c("time", "status"),
                   nName="cg")
        xDT[, "cg" := as.factor(cg)]
        xDT[, "strat" := xDTstr[, factor(strat)]]
        ## columns which are not strata
        xDTnotS1 <- seq.int(names(xDT))[!(grepl("^strat", names(xDT)))]
        res1 <- lapply(xDT[, seq.int(levels(strat))],
                       function(i) {
                           data.table::copy(xDT[as.integer(strat)==i, .SD, .SDcols=xDTnotS1])})
        ## drop unused levels
        res1 <- lapply(res1,
                       function(i) {
                           i[, "cg" := factor(cg)]})
        res1 <- lapply(res1,
                       ten,
                       abbNames=abbNames)
        ln1 <- data.table::data.table(
            "id" = xDT[, seq.int(levels(strat))],
            "longName" = xDT[, levels(strat)])
        data.table::setattr(res1,
                            name="longNames",
                            value=ln1)
        if (abbNames) {
            names(res1) <- ln1[, id]
        } else {
            names(res1) <- xDT[, levels(strat)]
        }
        data.table::setattr(res1,
                            name="class",
                            value=c("stratTen", class(res1)))
        data.table::setattr(res1, name="call", value=call)
        data.table::setattr(res1, name="abbNames", value=abbNames)
        return(res1)
    }
    if (!any(xNS1)) {
        mm1 <- data.table::copy(xDT)
        collapseDT(xDT,
                   except=c("time", "status"),
                   nName="cg")
        ten(x=xDT,
            abbNames=abbNames,
            mm=mm1,
            call=call)
    }
}
#' @rdname ten
#' @method ten data.table
#' @aliases ten.data.table
#' @export
#' 
ten.data.table <- function(x, ...,
                           abbNames=TRUE,
                           mm=NULL,
                           call=NULL){
    stopifnot(all(names(x) %in% c("time", "status", "cg")))
    partMatch(env1=environment(), ...)
    data.table::setkey(x, time, cg)
    ## number at risk
    x[, "n" := rev(seq.int(nrow(x)))]
    ## number at risk per covariate group
    x[, "ncg" := rev(seq.int(length(n))), by=cg]
    ## drop unused levels
    x[, "cg" := as.factor(cg)[, drop=TRUE]]
    ## long names
    ln1 <- data.table::data.table(
        "id" = x[, seq_along(levels(cg))],
        "longName" = x[, levels(cg)])
    if(abbNames) x[, "cg" := as.integer(cg)]
    data.table::setnames(x, c("t", colnames(x)[-1]))
    x[, "e" := sum(status), by=list(t, cg)]
    x[, "ncg" := max(ncg), by=list(t, cg)]
    x[, "n" := max(n), by=list(t)]
    x[, "status" := NULL]
    x <- x[!(duplicated(x)), ]
    data.table::setcolorder(x,
                            c("t", "e", "n", "cg", "ncg"))
    data.table::setkey(x, "cg")
    data.table::setattr(x, "class", c("ten", class(x)))
    setAttr(x,
            "shape"="long",
            "abbNames"=abbNames,
            "longNames"=ln1,
            "ncg"=nrow(ln1),
            "call"=call,
            "mm"=mm)
    return(x)
}
#' 
#' @rdname ten
#' @method ten ten
#' @aliases ten.ten
#' @export
#' 
ten.ten <- function(x, ...,
                    abbNames=NULL,
                    call=NULL){
    partMatch(env1=environment(), ...)
    if (attr(x, "shape")=="long") {
        return(asWide(x))
    } else {
        return(asLong(x))
    }
}
### helper functions
##
## partial matching with an ellipsis
## from environment env1
partMatch <- function(env1=NULL, ...){
    stopifnot(is.environment(env1))
    l1 <- as.list(substitute(list(...)))[-1L]
    n1 <- c("sh", "abb", "con")
    s1 <- sapply(n1, pmatch, names(l1))
    n2 <- c("shape", "abbNames", "contrasts.arg")
    names(s1) <- n2
    s1 <- s1[!is.na(s1)]
    for (i in seq_along(s1)){
        names(l1)[s1[i]] <- names(s1[i])
    }
    l1 <- l1[names(l1) %in% n2]
    for(i in seq_along(l1)){
        ## this isn't v. pretty...
        if (is.character(l1[[i]])){
            p1 <- paste0("env1$", names(l1)[i], " <- \"", l1[[i]], "\"")
        } else { 
            p1 <- paste0("env1$", names(l1)[i], " <- ", l1[[i]])
        }
        eval(parse(text=p1))
    }
}
## collapse/ paste a data table
## x = data.table
## except = columns to remain unmodified
## nName = new name for collapsed column
## returns the modified data.table
collapseDT <- function(x,
                       except=c("time", "status"),
                       nName="cg"){
    stopifnot(inherits(x, "data.table"))
    if (ncol(x)==1) {
        data.table::setnames(x, nName)
        return(invisible())
    }
    ## names in 'except'?
    toCollapse1 <- names(x)[!names(x) %in% except]
    x[, (nName) := paste(toCollapse1,
                         .SD,
                         sep="=",
                         collapse=", "),
      .SDcols=toCollapse1,
      by=seq.int(nrow(x))]
    toRemove1 <- which(names(x) %in% toCollapse1)
    if (length(toRemove1)) {
        data.table::set(x, j=toRemove1, value=NULL)
    }
    return(invisible())
}
## set attributes for a ten object (a data.table)
setAttr <- function(x, ...) UseMethod("setAttr")
setAttr.ten <- function(x, ...,
                        shape=NULL,
                        abbNames=NULL,
                        longNames=NULL,
                        ncg=NULL,
                        call=NULL,
                        mm=NULL){
    stopifnot(inherits(x, "ten"))
    ## can't use .Internal in a package...
    ## l1 <- .Internal(ls(envir=environment(), all.names=TRUE))
    l1 <- ls()
    l1 <- l1[!grepl("x", l1)]
    for(i in seq_along(l1)){
        data.table::setattr(x,
                            name=l1[i],
                            value=eval(as.name(l1[i])))
    }
    return(x)
}
## for R CMD check
n <- status <- strat <- time <- NULL
