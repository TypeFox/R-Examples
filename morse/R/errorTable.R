errorTableCreate <- function() {
  t <- data.frame(id = character(0), msg = character(0), stringsAsFactors = FALSE)
  class(t) <- c("errorTable", "data.frame")
  t
}

errorTableAppend <- function(...) {
  u <- rbind(...)
  class(u) <- c("errorTable", "data.frame")
  u
}

errorTableAdd <- function(t, id, msg) {
  newlines <- data.frame(id = id, msg = msg, stringsAsFactors = FALSE)
  errorTableAppend(t, newlines)
}

errorTableSingleton <- function(id, msg) {
  errorTableAdd(errorTableCreate(), id, msg)
}

errorTableIsEmpty <- function(x)
  dim(x)[1] == 0
  
#' @export
print.errorTable <- function(x, ...) {
  if (errorTableIsEmpty(x)) {
    cat("No error detected.\n")
  }
  else {
    cat("Error(s):\n")
    for (m in x$msg) {
      cat(paste("\t",m,"\n",sep=""))
    }
  }
}
