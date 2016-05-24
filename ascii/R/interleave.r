# from http://tolstoy.newcastle.edu.au/R/help/06/03/22717.html

##' interleave
##'
##' @export
##' @keywords internal
##' @param ... ...
interleave <- function(...) {
  args <- list(...)
  args <- args[!sapply(args, is.null)]
  nargs <- length(args)

  ord <- list()
  for (i in 1:nargs) {
    ord <- c(ord, nargs*(1:length(args[[i]])) - nargs + i)
  }
  do.call("c", args)[order(unlist(ord))]
}
## interleave(rep(1, 5),rep(3, 8)) 
## interleave(1:4, 5:8)
## interleave(1:4, 5:8, 9:12)
## interleave(list(1, 2, 3, 4), list(5, 6, 7, 8))


##' as.list.matrix
##' 
##' @keywords internal
##' @param x x
##' @param byrow byrow
##' @param ... ...
as.list.matrix <- function(x, byrow = TRUE, ...) {
  margin <- 2
  if (byrow)
    margin <- 1

  lapply(apply(x, margin, list), function(x) x[[1]])
}


##' interleave.matrix
##'
##' @export
##' @keywords internal
##' @param ... ...
##' @param byrow byrow
interleave.matrix <- function(..., byrow = TRUE) {
  args <- list(...)
  args <- args[!sapply(args, is.null)]
  
  lists <- lapply(args, function(x) {
    as.list.matrix(x, byrow = byrow)
  })
  interlists <- do.call("interleave", lists)
  if (byrow)
    do.call("rbind", interlists)
  else
    do.call("cbind", interlists)
}

##' interleave.data.frame
##'
##' @export
##' @keywords internal
##' @param ... ...
##' @param byrow byrow
##' @param pretty.rownames pretty.rownames
interleave.data.frame <- function(..., byrow = TRUE, pretty.rownames = TRUE) {
  args <- list(...)
  args <- args[!sapply(args, is.null)]
  
  names_df <- lapply(args, names)
  if (byrow) {
    inter_names <- names(args[[1]])
    class_df <- lapply(args[[1]], class)
    if (pretty.rownames) {
      args_names <- names(args)
      real_args_names <- as.character(as.list(substitute(list(...)))[-1])
      if (is.null(args_names))
        args_names <- real_args_names
      args_names[args_names == ""] <- real_args_names[args_names == ""]
      for (i in 1:length(args)) {
        row.names(args[[i]]) <- paste(args_names[i], row.names(args[[i]]), sep = ": ")
      }
    }
  } else {
    inter_names <- unlist(do.call("interleave", names_df))
  }
  
  list_mat <- lapply(args, as.matrix)
  names(list_mat) <- NULL
  results <- suppressWarnings(data.frame(do.call("interleave.matrix", c(list_mat, byrow = byrow))))
  names(results) <- inter_names
  if (byrow) {
    for (i in 1:ncol(results)) {
      class(results[, i]) <- class_df[[i]]
    }
  }
  results
}
