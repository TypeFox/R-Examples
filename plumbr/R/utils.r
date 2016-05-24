## utils

"%||%" <- function(a, b) if(is.null(a)) b else a

.irreversible <- function(reason)
  stop("Reversal impossible due to ", reason, ".")
