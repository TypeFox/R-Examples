is.likert <- function(x) inherits(x, "likert")

as.likert <- function(x, ...) {
  if (is.likert(x)) return(x)
  UseMethod("as.likert")
}

as.likert.data.frame <- function(x, ...) {
  as.likert(data.matrix(x), ...)
}

as.likert.formula <- function(x, ...) {
  data <- list(...)$data
  stop("use one of the idioms
    as.likert(tapply( <appropriate arguments> ))
    as.likert(table( <appropriate arguments> ))",
       call.=FALSE)
  ## data(apple)
  ## as.likert(tapply(apple$pre, apple[,1:2], c))
  ## as.likert(table(apple[,1:2]))
  ##
  ## eventually I will need to automate this:
  ##    data[,deparse(x[[2]])]), data[, x[[3]]]
  ## where you have to manually unpack x[[3]].
}


as.likert.ftable <- function(x, ...) {
  as.likert(as.table(x), ...)
}

as.likert.table <- function(x, ...) {
  as.likert.matrix(x, ...)
}

as.likert.matrix <- function(x,
                             ReferenceZero=NULL,
                             ...,
                             rowlabel=NULL, collabel=NULL,
                             xlimEqualLeftRight=FALSE,
                             xTickLabelsPositive=TRUE,
                             padding=FALSE,
                             reverse.left=TRUE) {
  ## All the as.likert calls end here.  This one accepts and ignores '...'.
  ## All the others may use the arguments here in their ... arguments.
  if (any(x < 0))
    stop("Argument to the likert() function must be non-negative.",
         call.=FALSE)
  nc <- ncol(x)
  colorset <- ColorSet(nc, ReferenceZero)

  ndnx <- names(dimnames(x))
  if (is.null(dimnames(x))) dimnames(x) <- list(1:nrow(x), NULL)
  if (is.null(dimnames(x)[[1]])) dimnames(x)[[1]] <- 1:nrow(x)

  levels <- dimnames(x)[[2]]
  if (is.null(levels)) {
    levels <- LETTERS[1:nc]
    dimnames(x)[[2]] <- levels
  }

  if (!is.null(rowlabel)) names(dimnames(x))[1] <- rowlabel
  if (!is.null(collabel)) names(dimnames(x))[2] <- collabel

  ind.neg <- (1:nc)[colorset < 0]
  ind.pos <- (1:nc)[colorset > 0]
  ind.zero <- (1:nc)[colorset == 0]

  x.neg <- x[,ind.neg, drop=FALSE]
  x.pos <- x[,ind.pos, drop=FALSE]
  x.zero <- x[,ind.zero, drop=FALSE]
  x.zero.pos <- x.zero
  if(0 %in% colorset)
    dimnames(x.zero.pos)[[2]] <- paste(dimnames(x)[[2]][ind.zero], "Positive")

  positive.order <- order(rowSums(cbind(x.zero.pos/2, x.pos)))

  if (padding) { ## for mosaic
    x.neg.pad <- {tmp <- as.matrix(rowSums(x.neg))
                  if (ncol(x.zero) > 0) tmp <- tmp + x.zero/2
                  max(tmp) - tmp}
    dimnames(x.neg.pad)[[2]] <- "left padding"
    x.pos.pad <- {tmp <- as.matrix(rowSums(x.pos))
                  if (ncol(x.zero.pos) > 0) tmp <- tmp + x.zero.pos/2
                  max(tmp) - tmp}
    dimnames(x.pos.pad)[[2]] <- "right padding"
    ind.pad <- 1
  } else { ## for barchart
    x.neg.pad <- x[, 0, drop=FALSE]
    x.pos.pad <- x.neg.pad
    ind.pad <- 0
  }

  if (reverse.left) { ## for barchart
    x <- cbind(-x.zero/2,    -x.neg[, rev(seq_len(ncol(x.neg))), drop=FALSE], -x.neg.pad,
               x.zero.pos/2,  x.pos,                                           x.pos.pad)

    attr(x, "color.seq") <- c(ind.zero, rev(ind.neg), 0[ind.pad],
                              ind.zero,     ind.pos,  0[ind.pad])
  } else { ## for mosaic
    x <- cbind(x.neg.pad,    x.neg, x.zero/2,
               x.zero.pos/2, x.pos, x.pos.pad)

    attr(x, "color.seq") <- c(0[ind.pad], ind.neg, ind.zero,
                              ind.zero, ind.pos, 0[ind.pad])
  }

  attr(x, "positive.order") <- positive.order

  names(dimnames(x)) <- ndnx
  attr(x, "nlevels") <- nc
  attr(x, "original.levels") <- levels
  attr(x, "xlimEqualLeftRight") <- xlimEqualLeftRight
  attr(x, "xTickLabelsPositive") <- xTickLabelsPositive
  attr(x, "ReferenceZero") <- ReferenceZero
  class(x) <- c("likert", class(x))
  x
}


## environment(as.likert.matrix) <- environment(plot.likert)
## assignInNamespace(x, value, ns, pos = -1,
##                   envir = as.environment(pos))
## assignInNamespace("as.likert.matrix", as.likert.matrix, "HH")


as.likert.default <- function(x, ...) {
  ## simple vector because anything else got dispatched elsewhere
  x <- t(x)
  if (is.null(dimnames(x))) dimnames(x) <- list("", 1:length(x))
  if (is.null(dimnames(x)[[1]])) dimnames(x)[[1]] <- ""
  as.likert(x, ...)
}

is.likertCapable <- function(x, ...) {
  is.numeric(x) ||
  is.table(x) ||
  inherits(x, "ftable") ||
  ("package:vcd" %in% search() && is.structable(x)) ||
  is.data.frame(x) ||
  is.listOfNamedMatrices(x)
}


as.likert.listOfNamedMatrices <- function(x, ...) {
  result <- sapply(x, as.likert, simplify=FALSE, ...)
  class(result) <- c("likert", "list")
  result
}

as.likert.array <- function(x, ...)
  as.likert(as.listOfNamedMatrices(x), ...)

rev.likert <- function(x) {
  ## Reverses the rows of a matrix "likert" object,
  ## or each item within a list "likert" object,
  ## and retains all attributes.

  if (is.matrix(x)) {
    z <- x
    z[] <- if (nrow(x)) x[nrow(x):1L, , drop=FALSE] else x
    dimnames(z)[[1]] <- rev(dimnames(x)[[1]])
    z
  }
  else {
    sapply(x, rev, simplify=FALSE)
  }
}


## This simplified function appears in the JSS submitted paper.
## It is here for discussion but it not exported.
## The user can access it with HH:::as.likert.simplified.odd
as.likert.simplified.odd <- function(x,                          ##  1
       nc=ncol(x), colorset=(1:nc)-(nc+1)/2) {                   ##  2
    ind.neg <- rev((1:nc)[colorset < 0])                         ##  3
    ind.pos <- (1:nc)[colorset > 0]                              ##  4
    ind.zero <- (1:nc)[colorset == 0]                            ##  5
    x <- cbind(-x[, ind.zero, drop = FALSE]/2,                   ##  6
               -x[, ind.neg, drop = FALSE],                      ##  7
                x[, ind.zero, drop = FALSE]/2,                   ##  8
                x[, ind.pos, drop = FALSE])                      ##  9
    attr(x, "color.seq") <- c(ind.zero, ind.neg, ind.pos)        ## 10
    pos.columns <- seq(to = ncol(x),                             ## 10
         length = length(c(ind.zero, ind.pos)))                  ## 12
    attr(x, "positive.order") <- order(apply(x[, pos.columns,    ## 13
        drop = FALSE], 1, sum))                                  ## 14
    x                                                            ## 15
}                                                                ## 16

## source("c:/HOME/rmh/HH-R.package/HH/R/as.likert.R")
