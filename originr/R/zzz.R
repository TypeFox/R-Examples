mssg <- function(v, ...) if (v) message(...)

orc <- function(l) Filter(Negate(is.null), l)
