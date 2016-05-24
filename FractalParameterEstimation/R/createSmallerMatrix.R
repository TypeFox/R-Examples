createSmallerMatrix <-
function(givenMatrix) {
  if(missing(givenMatrix)) 
    stop("Missing data")
  ## calculates ramification of given matrix
  ## deducts 1 for new ramification
  ramification = calcRamification(dim(givenMatrix)[[1]]) - 1
  ## create new matrix
  dimension = potence(3, ramification)*potence(3, ramification)
  x = matrix(as.integer(rep(99,dimension)), potence(3, ramification), potence(3, ramification))
  return(x)
}
