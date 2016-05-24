as.character.GEPhenotype <- function(x, ..., simplify=TRUE) {
  expr = x$parsed
    
  if (simplify) {
    expr = parse(text=as.character(expr))
  }

  return (as.character(expr))
}

as.expression.GEPhenotype <- function(x, ...) {
  if (x$type == "NT") {
    warning("Non-Terminal Sequence")
    NULL
  } else {
    x$parsed
  }
}

print.GEPhenotype <- function(x, ..., simplify=TRUE) {
  if (x$type == "NT") {
    cat("Non-Terminal Sequence:\n", x$expr, "\n", ...)
    #print(NA, ...)
  } else {
    cat(as.character(x, simplify=simplify), "\n", ...)
  }
}
