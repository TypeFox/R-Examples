#' Extended Tabular Operations for Both matrix and big.matrix Objects
#'
#' @rdname bigtabulate
#' @useDynLib bigtabulate
#' @aliases bigtabulate bigsplit bigtable bigtsummary
#' @description This package extends the \pkg{bigmemory} package, but the 
#' functions may also be used with traditional \R \code{matrix} and 
#' \code{data.frame} objects. The function \code{\link{bigtabulate}} is 
#' exposed, but we expect most users will prefer the higher-level functions 
#' \code{bigtable}, \code{bigtsummary}, and \code{bigsplit}. Each of these
#' functions provides functionality based on a specified conditional 
#' structure.  In other words, for every cell of a (possibly multidimensional) 
#' contingency table, they provide (or tabulate) some useful conditional 
#' behavior (or statistic(s)) of interest.  At the most basic level, this 
#' provides an extremely fast and memory-efficient alternative to 
#' \code{\link{table}} for matrices and data frames.
#' @param x a \code{\link[bigmemory]{big.matrix}} or a 
#' \code{\link{data.frame}} or a \code{\link{matrix}}.
#' @param ccols a vector of column indices or names specifying which 
#' columns should be used for conditioning (e.g. for building a contingency 
#' table or structure for tabulation).
#' @param breaks a vector or list of \code{length(ccols)}.  If a vector,
#' \code{NA} indicates that the associated column should be treated like a
#' factor (categorical variable), while an integer value indicates that the
#' range of the associated column should be broken into a specified number of
#' evenly-spaced bins (histogram-like).  If a list, \code{NA} triggers the
#' factor-like handling, a single number triggers bin-like behavior, while a
#' triplet (min,max,breaks) indicates that the bin-like behavior should be on
#' a restricted range rather than on the range of data for that column.  See
#' \code{\link[biganalytics]{binit}} for similar specification of this option.
#' @param table if \code{TRUE}, a list of table counts will be returned.
#' @param useNA whether to include extra '\code{NA}' levels in the table.
#' @param summary.cols column(s) for which table summaries will be calculated.
#' @param summary.na.rm if \code{TRUE}, \code{NA}s are removed from table
#' summary calculations.
#' @param splitcol if \code{NA}, the indices which correspond to
#' table-levels are returned.  If numeric, the corresponding column
#' values will be returned in a list corresponding to table-levels.  If
#' \code{NULL}, then there is no splitting at all.
#' @param splitret if \code{"list"}, the \code{splitcol} value is returned
#' as a list.  When \code{splitcol} is \code{NA}, \code{splitret} may
#' be \code{"vector"}.  Finally, \code{"sparselist"} may be a useful option
#' when the full-blown splitting structure has many unrepresented "cells";
#' this is like using the \code{drop=TRUE} option to \code{\link{split}}.
#' @param cols with \code{bigtsummary}, which column(s) should be conditionally
#' summarized?  This (or these) will be passed on as \code{summary.cols}.
#' @param na.rm an obvious option for summaries.
#' @return array-like object(s), each similar to what is returned by
#' \code{\link{tapply}} and the associated \R functions.
#' @details This package concentrates on conditional stuctures and calculations,
#' much like \code{\link{table}}, \code{\link{tapply}}, and \code{\link{split}}.
#' The functions are juiced-up versions of the base \R functions;
#' they work on both regular \R matrices and data frames, but are specialized
#' for use with \pkg{bigmemory} and (for more advanced usage) \pkg{foreach}.
#' They are particularly fast and memory-efficient.  We have found that
#' \code{bigsplit} followed by \code{\link{lapply}} or \code{\link{sapply}}
#' can be particularly effective, when the subsets produced by the split
#' are of reasonable size.  For intensive calculations, subsequent use of
#' \code{foreach} can be helpful (think: parallel apply-like behavior).
#' 
#' When \code{x} is a \code{matrix} or a \code{data.frame}, some additional
#' work may be required.  For example, a character column of a \code{data.frame}
#' will be converted to a \code{\link{factor}} and then coerced to numeric
#' values (factor level numberings).
#' 
#' The conditional structure is specified via \code{ccols} and \code{breaks}.
#' This differs from the design of the base \R functions but is at the root
#' of the gains in speed and memory-efficiency.  The \code{breaks} may seem
#' distracting, as most users will simply condition on categorical-like columns.
#' However, it provides the flexibility to \dQuote{bin} \dQuote{continuous},
#' column(s) much like a histogram.  See \code{\link[biganalytics]{binit}} for
#' another example
#' of this type of option, which can be particularly valuable with massive 
#' data sets.
#' 
#' A word of caution: if a \dQuote{continuous} variable is not \dQuote{binned},
#' it will be treated like a factor and the resulting conditional structure will
#' be large (perhaps immensely so).
#' The function uses left-closed intervals [a,b) for the "binning" behavior,
#' when specified, except in the right-most bin, where the interval is entirely
#' closed.
#' 
#' Finally, \code{bigsplit} is somewhat more general than \code{split}.
#' The default behavior (\code{splitcol=NA})
#' returns a split of \code{1:nrow(x)} as a list
#' based on the specified conditional structure.  However, it may also
#' return a vector of cell (or category) numbers.  And of course it may
#' conduct a split of \code{x[,splitcol]}.
#' @export
#' @examples
#' data(iris)
#' 
#' # First, break up column 2 into 5 groups, and leave column 5 as a
#' # factor (which it is).  Note that iris is a data.frame, which is
#' # fine.  A matrix would also be fine.  A big.matrix would also be fine!
#' bigtable(iris, ccols=c(2, 5), breaks=list(5, NA))
#' 
#' iris[,2] <- round(iris[,2]) # So columns 2 and 5 will be factor-like
#'                             # for convenience in these examples, below:
#' 
#' ans1 <- bigtable(iris, c(2, 5))
#' ans1
#' # Same answer, but with nice factor labels from table(), because
#' # table() handles factors.  bigtable() uses the numeric factor
#' # levels only.
#' table(iris[,2], iris[,5])
#' 
#' # Here, our formulation is simpler than split's, and is faster and
#' # more memory-efficient:
#' ans2 <- bigsplit(iris, c(2, 5), splitcol=1)
#' ans2[1:3]
#' split(iris[,1], list(col2=factor(iris[,2]), col5=iris[,5]))[1:3]
#' @import bigmemory biganalytics
bigtabulate <- function(x,
                        ccols, breaks=vector("list", length=length(ccols)),
                        table=TRUE, useNA="no",
                        summary.cols=NULL, summary.na.rm=FALSE,
                        splitcol=NULL, splitret="list") {

  if (!is.matrix(x) && !is.data.frame(x)) {
    if (class(x)!="big.matrix")
      stop("bigtabulate requires matrix, data.frame, or big.matrix objects.")
  }

  if (is.data.frame(x)) {
    for (i in 1:ncol(x)) {
      if (is.character(x[,i])) x[,i] <- factor(x[,i])
      if (is.factor(x[,i])) x[,i] <- as.integer(x[,i])
    }
    x <- as.matrix(x)
  }

  # Check and prepare ccols
  if (is.logical(ccols)) ccols <- which(ccols)
  if (length(ccols)!=length(breaks))
    stop("length(ccols) must equal length(breaks).")
  if (!is.numeric(ccols) & !is.character(ccols))
    stop("column indices must be numeric or character vectors.")
  if (is.character(ccols))
    if (is.null(colnames(x))) stop("column names do not exist.")
    else ccols <- bigmemory:::mmap(ccols, colnames(x))

  # Prepare breaks: could be a vector of length(ccols) of numbers of
  # breaks (or NA), assumed to span the ranges of the variables; or a list of
  # the same length containing a mixture of numbers of breaks (or NA) or triplets
  # of (min, max, breaks).  The result of the preparation is a matrix
  # with 3 rows and length(ccols) columns of (min, max, breaks) values (possibly NA).
  # NA indicates factor-like handling of the variable for the tabulations.
  breaks[sapply(breaks, is.null)] <- NA
  breakm <- matrix(NA, 3, length(breaks))
  if (is.numeric(breaks)) {
    if (is.matrix(x)) {
      breakm[1,!is.na(breaks)] <- apply(x[,ccols[!is.na(breaks)], drop=FALSE], 2, min, na.rm=TRUE)
      breakm[2,!is.na(breaks)] <- apply(x[,ccols[!is.na(breaks)], drop=FALSE], 2, max, na.rm=TRUE)
    } else {
      breakm[1,!is.na(breaks)] <- biganalytics::colmin(x, ccols[!is.na(breaks)], na.rm=TRUE)
      breakm[2,!is.na(breaks)] <- biganalytics::colmax(x, ccols[!is.na(breaks)], na.rm=TRUE)
    }
    breakm[3,] <- breaks
  }
  if (is.list(breaks)) {
    for (i in which(!sapply(breaks, is.na))) {
      if (length(breaks[[i]])==1) {
        if (is.matrix(x)) { 
          breakm[1,i] <- min(x[,ccols[i]], na.rm=TRUE)  
          breakm[2,i] <- max(x[,ccols[i]], na.rm=TRUE)
        } else {
          breakm[1,i] <- biganalytics::colmin(x, ccols[i], na.rm=TRUE)
          breakm[2,i] <- biganalytics::colmax(x, ccols[i], na.rm=TRUE)
        }
        breakm[3,i] <- breaks[[i]]
      } else {
        breakm[,i] <- breaks[[i]]
      }
    }
  }

  if (!is.logical(table)) stop("table must be logical.")
  table.useNA <- -1
  if (useNA=="no") table.useNA <- as.integer(0)
  if (useNA=="ifany") table.useNA <- as.integer(1)
  if (useNA=="always") table.useNA <- as.integer(2)
  if (table.useNA==-1) stop("invalid argument to useNA.")

  if (!is.logical(summary.na.rm)) stop("summary.na.rm must be logical.")
  if (is.logical(summary.cols)) summary.cols <- which(summary.cols)
  if (!is.numeric(summary.cols) && !is.character(summary.cols) &&
    !is.null(summary.cols)) {
    stop(paste("summary column indices must be numeric, logical,",
       "or character vectors."))
  }
  if (is.character(summary.cols))
    if (is.null(colnames(x))) stop("column names do not exist.")
    else summary.cols <- bigmemory:::mmap(summary.cols, colnames(x))
  if (!is.null(splitcol)) {
    if (!is.na(splitcol)) {
      if (is.logical(splitcol)) splitcol <- which(splitcol)
      if (!is.numeric(splitcol) & !is.character(splitcol))
        stop("splitcol must be numeric, logical, or character specifying one column, or NA or NULL.")
      if (is.character(splitcol))
        if (is.null(colnames(x))) stop("column names do not exist.")
        else splitcol <- bigmemory:::mmap(splitcol, colnames(x))
      if (length(splitcol)!=1) stop("splitcol must identify a single column or be NA or NULL.")
      splitcol <- as.numeric(splitcol)
    }
  }

  if (splitret!="vector" && splitret!="list" && splitret!="sparselist")
    stop("splitret must be 'vector' or 'list' or 'sparselist'")
  splitlist <- 0 # Was FALSE; this indicates vector return possibility
  if (splitret=="list") splitlist <- 1 # Was TRUE
  if (splitret=="sparselist") splitlist <- 2 # New option
  if (splitret=="vector" && is.numeric(splitcol))
    stop("splitting a specified column must return a list")
  splitlist <- as.integer(splitlist)

  # splitcol=NULL	Don't return any map type of anything.
  # splitcol=NA		Essentially split 1:nrow(x)
  # splitcol=a column   Split this single column.
  # splitlist=2: an option for the sparse list representation.
  # splitlist=1 by default: the return is a list of either split 1:nrow(x) or col entries
  # splitlist=0: only valid if splitcol==NA, in which case a vector of as.numeric(factor) entries
  summary <- is.numeric(summary.cols)
  if (is.numeric(splitcol) && splitlist==0)
    stop("vector split structure is not allowed on a column")

  if (!is.matrix(x)) {
    ans <- .Call("BigMatrixTAPPLY", x@address, as.numeric(ccols), 
                 as.numeric(breakm),
                 table, table.useNA,
                 summary, as.numeric(summary.cols), 
                 summary.na.rm, splitcol, splitlist,
                 PACKAGE="bigtabulate")
  } else {
    if (is.integer(x)) {
      ans <- .Call("RIntTAPPLY", x, as.numeric(ccols), as.numeric(breakm),
                   table, table.useNA,
                   summary, as.numeric(summary.cols), 
                   summary.na.rm, splitcol, splitlist,
                   PACKAGE="bigtabulate")
    } else {
      ans <- .Call("RNumericTAPPLY", x, as.numeric(ccols), as.numeric(breakm),
                   table, table.useNA,
                   summary, as.numeric(summary.cols), 
                   summary.na.rm, splitcol, splitlist,
                   PACKAGE="bigtabulate")
    }
  }

  # The return will always contain
  # - ans$levels, a list of length(ccols) of factor levels possibly plus "NA"
  #
  # It will contain at least one of the following:
  # - ans$table:	vector of length prod(dim())
  # - ans$summary:	list of length prod(dim()) of cell summary matrices with 5 columns
  # - ans$split:	list of length prod(dim()) containing the split or map result;
  #                     nothing returned if is.null(splitcol), so don't return anything from C++
  #                     in that case.  Or a vector of the factor levels.

  dn <- lapply(ans$levels, function(x) { x[is.na(x)] <- "NA"; return(x) })
  ans$levels <- NULL
  if (table) ans$table <- array(ans$table, dim=sapply(dn, length), dimnames=dn)
  if (summary){
     ans$summary <- array(ans$summary, dim=sapply(dn, length), dimnames=dn)
  }
  #if (!is.null(splitcol)) {
  #  names(ans$split) <- unlist(dn) # This currently only works for 1-factor splits.
  #}

  if (length(ans)==1) return(ans[[1]])
  return(ans)

}

#' @export
#' @rdname bigtabulate
bigsplit <- function(x, ccols,
                     breaks=vector("list", length=length(ccols)), useNA="no", 
                     splitcol=NA, splitret="list") {

  return(bigtabulate(x, ccols=ccols, breaks=breaks,
                     table=FALSE, useNA=useNA,
                     splitcol=splitcol, splitret=splitret))
}

#' @export
#' @rdname bigtabulate
bigtable <- function(x, ccols,
                     breaks=vector("list", length=length(ccols)),
                     useNA="no") {

  return(bigtabulate(x, ccols=ccols, breaks=breaks,
                     table=TRUE, useNA=useNA,
                     splitcol=NULL))
}

#' @export
#' @rdname bigtabulate
bigtsummary <- function(x, ccols,
                        breaks=vector("list", length=length(ccols)), useNA="no",
                        cols, na.rm=FALSE) {

  return(bigtabulate(x, ccols=ccols, breaks=breaks,
                     table=FALSE, useNA=useNA,
                     summary.cols=cols, summary.na.rm=na.rm))

}

#
# April 25, 2010: we decided not to include bigaggregate() at this point.
# It just wasn't adding that much, and the performance is poor for all but
# the largest examples.
#
#bigaggregate <- function(x, stats, usesplit=NULL,
#                         ccols=NA, breaks=vector("list", length=length(ccols)), 
#                         useNA="no", distributed=FALSE, rettype="celllist", 
#                         simplify=TRUE) {
#  if (is.null(usesplit)) {
#    usesplit <- bigsplit(x, ccols=ccols, breaks=breaks, useNA=useNA, 
#      splitcol=NA, splitret="list")
#  }
#
#  # At this point I have usesplit, which is the map.  Everything else is much like I had
#  # previously in commented code, below.
#
#  if (is.data.frame(x)) {
#    for (i in 1:ncol(x)) {
#      if (is.character(x[,i])) x[,i] <- factor(x[,i])
#      if (is.factor(x[,i])) x[,i] <- as.integer(x[,i])
#    }
#    x <- as.matrix(x)
#  }
#
#  require(foreach)
#  if (is.null(getDoParName())) {
#    registerDoSEQ() # A little hack to avoid the foreach warning 1st time.
#  }
#  if (!is.list(stats[[1]])) stats <- list(stats=stats)
#  if (is.null(names(stats))) stop("stats must be a named list")
#  if (length(unique(names(stats)))!=length(stats))
#    stop("names of stats list must be unique")
#
#  if (!is.null(getDoParName()) && getDoParName()!="doSEQ") {
#    require(bigmemory)
#    if (distributed) bf <- "" 
#    else bf <- NULL
#    if (is.matrix(x)) {
#      x <- as.big.matrix(x, backingfile=bf)
#      warning("Temporary shared big.matrix created for parallel calculations.")
#    }
#    if (!is.shared(x) && !is.filebacked(x)) {
#      x <- deepcopy(x, backingfile=bf)
#      warning("Temporary shared big.matrix created for parallel calculations.")
#    }
#  }
#
#  # Now prepare the arguments.
#  for (i in 1:length(stats)) {
#    thisname <- names(stats)[i]
#    args <- stats[[i]]
#    if (!is.list(args))
#      stop(paste("stats element", thisname, "needs to be list."))
#    if (!is.function(args[[1]]) && !is.character(args[[1]])) {
#      stop(paste("first argument of stats element", thisname, 
#        "needs to be a function."))
#    }
#    if (is.character(args[[2]])) {
#      if (is.null(colnames(x))) stop("column names do not exist.")
#      else args[[2]] <- bigmemory:::mmap(args[[2]], colnames(x))
#    }
#    if (!is.numeric(args[[2]])) args[[2]] <- as.numeric(args[[2]])
#    stats[[i]] <- args
#    names(stats)[i] <- thisname
#  }
#
#  # Here, process the chunks of data.
#  xdesc <- if (!is.matrix(x)) describe(x) else NULL
#  fans <- foreach(i=usesplit) %dopar% {
#    if (is.null(i)) {
#      temp <- as.list(rep(NA, length(stats)))
#      names(temp) <- names(stats)
#      return(temp)
#    }
#    if (!is.null(xdesc)) x <- attach.big.matrix(xdesc)
#    temp <- vector("list", length=0)
#    for (j in names(stats)) {
#      farg <- stats[[j]]
#      tempname <- names(formals(farg[[1]]))[1]
#      if (is.character(farg[[1]])) farg[[1]] <- as.symbol(farg[[1]])
#      farg[[2]] <- x[i,farg[[2]],drop=FALSE]
#      if (!is.null(tempname)) names(farg)[2] <- tempname
#      else names(farg)[2] <- ""
#      mode(farg) <- "call"
#      temp[[j]] <- eval(farg)
#    }
#    rm(farg)
#    gc()
#    
#    return(temp)
#  }
#
#  #temp <- array(temp, dim=sapply(dn, length), dimnames=dn)
#
#  if (rettype=="statlist") {
#    # Provide list of length(stats) of arrayed answers.
#    z <- NULL
#    for (j in names(stats)) {
#      res <- lapply(temp, function(x) return(x[[j]]))
#      if (all(sapply(res, length)==1) && simplify)
#        res <- array(unlist(res), dim=sapply(dn, length), dimnames=dn)
#      else {
#        # Here, there could be some empty cells with single NA values that need replication:
#
##        if (length(unique(sapply(res, length)))==2) {
##          nc <- max(unique(sapply(temp, length)), na.rm=TRUE)
##          usenames <- names(temp[[which(sapply(temp, length)==nc)[1]]])
##          for (k in which(sapply(temp, length)==1)) {
##            if (is.na(temp[[k]])) temp[[k]] <- as.numeric(rep(NA, nc))
##            names(temp[[k]]) <- usenames
##          }
##        }
#
#        res <- array(temp, dim=sapply(dn, length), dimnames=dn)
#      }
#      z[[j]] <- temp
#    }
#  }
#
#  z[is.null(z)] <- NULL
#
#  return(fans)
#
#}

