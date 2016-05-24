#' @name lload
#' @author Sven E. Templer
#' @title Load RData Objects to a List
#' @description
#' Load multiple .RData files and return a (simplified) list.
#' @param path Character string with the path, as used in \link{list.files}.
#' @param pattern A regular expression for file name patterns, as used in
#' \link{list.files}.
#' @param recursive Logical. Search the path recursive.
#' @param simplify Logical, unlist when there are only unique object names.
#' @param verbose Logical. Print information on screen about loading process.
#' @return
#' Returns a list of length n, when there are n data files loaded. All objects
#' are stored in sublists. Names are according to files, and names of sublists
#' to objects per file. If simplified, the list is of length m, when there are
#' m objects in total loaded.
#' @seealso
#' \link{load}

#' @export lload
lload <- function (path = ".", pattern = ".RData", recursive = FALSE, simplify = TRUE, verbose = TRUE) {
  rds <- list.files(path, pattern, full.names = FALSE, recursive = recursive)
  lst.rds <- list()
  if (verbose) cat("Loading", length(rds), "data files ...")
  for (rd in rds) {
    e <- new.env()
    obj <- try(load(file.path(path, rd), e))
    if (class(obj) != "try-error") {
      lst.obj <- list()
      for (i in obj) {
        lst.obj[[i]] <- get(i, envir = e)
      }
      if (length(obj))
        lst.rds[[rd]] <- lst.obj
    }
  }
  if (verbose) cat("Done.\n")
  # unlist
  if (simplify) {
    if (verbose) cat("Simplifying ... ")
    objnames <- sapply(lst.rds, names)
    if (length(unique(objnames)) == length(objnames)) {
      names(lst.rds) <- NULL
      lst.rds <- unlist(lst.rds, recursive = FALSE)
      if (verbose) cat("Done.\n")
    } else {
      if (verbose) cat("Duplicates found ... Skipped.\n")
    }
  }
  return(lst.rds)
}