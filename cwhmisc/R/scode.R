scode <- function(p) {  # lifted from stats::printCoefmat
  p <- abs(p)
  if (p <= 0.001) return("***") else
    if (p <= 0.01) return("**") else
      if (p <= 0.05) return("*") else
        if (p <= 0.1) return(".") else
          return(" ")
}
