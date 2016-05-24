as.matrix.listOfNamedMatrices <- function(x, abbreviate=TRUE, minlength=4, ...) {
  result <- is.listOfNamedMatrices(x, xName=deparse(substitute(x)))
  if (!result) {
    stop(attr(result,"reason"))
  }
  ggg <- names(x)
  if (abbreviate) ggg <- abbreviate(ggg, minlength=minlength)

  nnn <- sapply(x,
                function(xi, abbreviate=abbreviate, minlength=minlength) {
                  nn <- rownames(xi)
                  if (abbreviate) nn <- abbreviate(nn, minlength=minlength)
                },
                abbreviate=abbreviate,
                minlength=minlength,
                simplify=FALSE)
  mmm <- do.call("rbind", x)
  rownames(mmm) <- paste(rep(ggg, sapply(x, nrow)),
                         unlist(nnn, use.names=FALSE),
                         sep=" ")
  attr(mmm, "Subtables.Rows") <- sapply(x, rownames, simplify=FALSE)
  mmm
}


is.listOfNamedMatrices <- function(x, xName=deparse(substitute(x))) {
  force(xName)
  result <- inherits(x, "listOfNamedMatrices")
  if (result) return(result)
  result <- is.list(x) && !is.data.frame(x)
  if (!result) {
    attr(result, "reason") <- paste(xName, "is not a list")
    return(result)
  }
  result <- !is.null(names(x))
  if (!result) {
    attr(result, "reason") <- paste("items in", xName, "are not named.")
    return(result)
  }
  for (nxi in names(x)) { ## convert vectors to single-row matrices
    xi <- x[[nxi]]
    if (is.numeric(xi) && (is.null(dim(xi)) || length(dim(x)) == 1)) x[[nxi]] <- t(xi)
  }
  result <- all(sapply(x, function(x) length(dim(x))) == 2)
  if (!result) {
    attr(result, "reason") <- paste("At least one item in", xName, "has more than two dimensions.")
    return(result)
  }
  for (nxi in names(x)) { ## verify that any data.frames have only numeric columns
    xi <- x[[nxi]]
    if (is.data.frame(xi))
      result <- (all(sapply(xi, is.numeric)))
    if (!result) {
      attr(result, "reason") <-
        paste("At least one item in", xName, "is a data.frame with a non-numeric column.")
      return(result)
    }
  }
  result <- all(sapply(x, ncol) == ncol(x[[1]]))
  if (!result) {
    attr(result, "reason") <- paste("Not all items in", xName, "have the same number of columns.")
    return(result)
  }
  result
}

as.listOfNamedMatrices <- function(x, xName=deparse(substitute(x)), ...)
  UseMethod("as.listOfNamedMatrices")


as.listOfNamedMatrices.list <- function(x, xName=deparse(substitute(x)), ...) {
  force(xName)
  result <- is.listOfNamedMatrices(x, xName=xName)
  if (!result) {
    stop(attr(result,"reason"))
  }
  if (!inherits(x, "listOfNamedMatrices"))
    class(x) <- c("listOfNamedMatrices", class(x))
  x
}

as.listOfNamedMatrices.matrix <- function(x, xName=deparse(substitute(x)), ...) {
  force(xName)
  tmp2 <- data.matrix(x)
  dim(tmp2) <- c(dim(tmp2)[1], 1, dim(tmp2)[2])
  dimnames(tmp2) <- list(dimnames(x)[[1]],
                         "nonsense",
                         dimnames(x)[[2]])
  tmp3 <- as.listOfNamedMatrices(aperm(tmp2, c(2,3,1)), xName=xName, ...)
  tmp4 <- sapply(names(tmp3),
                 function(x) {
                   dimnames(tmp3[[x]])[[1]] <- x
                   tmp3[[x]]},
                 simplify=FALSE)
  class(tmp4) <- class(tmp3)
  tmp4
}

as.listOfNamedMatrices.data.frame <- function(x, xName=deparse(substitute(x)), ...) {
  force(xName)
  as.listOfNamedMatrices(data.matrix(x), xName=xName, ...)
}

as.listOfNamedMatrices.MatrixList <- function(x, xName=deparse(substitute(x)), ...) {
  force(xName)
  NextMethod("as.listOfNamedMatrices")
}

as.listOfNamedMatrices.array <- function(x, xName=deparse(substitute(x)), ...) {
  force(xName)
  as.listOfNamedMatrices(as.MatrixList(x), xName=xName)
}


as.data.frame.listOfNamedMatrices <- function(x, ...) {
  xName <- deparse(substitute(x))
## old.warn <- options(warn=1)
  warning(paste("##", xName, "is a 'listOfNamedMatrices' and will not be converted to a data.frame."),
          call.=FALSE)
## recover()
## ## if (sys.nframe() > 10) { ## inside Rcmdr
## ##   Rcmdr::doItAndPrint(paste("##", xName, "remains a 'listOfNamedMatrices'.
## ## It's items are not variables in a data.frame.
## ## You may ignore the messages:
## ## in the Rcmdr Messages Window:
## ##   ERROR: the dataset", xName, "is no longer available.
## ## in the R Console:
## ##   Error in get(dataSet, envir = .GlobalEnv) : invalid first argument."))
## ## }
## options(old.warn)
  x
}

print.listOfNamedMatrices <- function(x, ...) {
  cat("'listOfNamedMatrices'.\n")
  print(as.matrix(x), ...)
  invisible(x)
}

`[.listOfNamedMatrices` <- function(x, ...) {
  result <- NextMethod("[")
  class(result) <- class(x)
  result
}


as.MatrixList <- function(x)
  UseMethod("as.MatrixList")

as.MatrixList.array <- function(x) {
  if (is.null(dimnames(x)) || any(sapply(dimnames(x), is.null))) stop("The object must have dimnames.")
  ldx <- length(dim(x))
  xa <- lapply(apply(x, 3:ldx, function(x) list(x)), `[[`, 1)
  dim(xa) <- dim(x)[-(1:2)]
  dimnames(xa) <- dimnames(x)[-(1:2)]
  names(dimnames(xa)) <- names(dimnames(x))[-(1:2)]
  if (is.null(names(xa))) { ## getting here means ldx > 3
    nxa <- outer(dimnames(x)[[3]], dimnames(x)[[4]], "paste", sep=".")
    if (ldx >= 5) {
      for (i in 5:ldx)
      nxa <- outer(nxa, dimnames(x)[[i]], "paste", sep=".")
    }
    names(xa) <- nxa
  }
  class(xa) <- c("MatrixList", "list", class(xa))
  xa
}
## environment(as.MatrixList.array) <- environment(plot.likert)

print.MatrixList <- function(x, ...) {
    cat("'MatrixList'.\n")
  print(as.listOfNamedMatrices(x), ...)
  invisible(x)
}


## source("c:/HOME/rmh/HH-R.package/HH/R/as.matrix.listOfMatrices.R")
