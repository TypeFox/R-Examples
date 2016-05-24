"is.inf" <- function(x) is.infinite(x)
"is.number" <- function(x) (is.numeric(x) || is.complex(x)) & !is.na(x)
"is.orderable" <- function(x) !is.na(x)