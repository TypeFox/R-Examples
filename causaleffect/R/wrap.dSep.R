wrap.dSep <- 
function(amat, first, second, cond) {
  if (identical(first, second)) return(FALSE)
  if (length(first) == 0 | length(second) == 0) {
    return(TRUE)
  }
  else return(dSep(amat, first, second, cond))
}