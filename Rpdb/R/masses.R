masses <- function(...)
  UseMethod("masses")

masses.default <- function(x, ...) {
  if(is.character(x)) return(elements$mass[match(x, elements$symb)])
  else if(is.numeric(x) & x == round(x) ) return(elements$mass[match(x, elements$num)])
  else stop("Bad argument: 'x' must be a character or an integer vector")
}

masses.pdb <- function(x, ...) {
  x <- toSymbols(x$atoms$elename)
  return(masses(x))
}