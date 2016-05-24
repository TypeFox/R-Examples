### 3 components - ternary
toTernary <- function(abc){
  sqrt3     <- 1.732050807568877293527446341505872366942805253810380628055806979

  return(cbind(
    x = (abc[,1L] + 2.0*abc[,3L]) / sqrt3,
    y = abc[,1L]
  ))

}

toTernaryVectors <- function(c1, c2, c3){ return(toTernary(cbind(c1, c2, c3))) }

### 4 components - quaternary
toQuaternary <- function(abcd){
  sqrt3     <- 1.732050807568877293527446341505872366942805253810380628055806979

  return(cbind(
    x = (abcd[,1L] + 2.0*abcd[,3L] + abcd[,4L]) / sqrt3,
    y = abcd[,1L] + abcd[,4L]/3.0,
    z = abcd[,4L]
  ))

}

toQuaternaryVectors <- function(c1, c2, c3, c4){ return(toQuaternary(cbind(c1, c2, c3, c4))) }

### "generic" function
toSimplex <- function(x){
  # checks
  if(is.null(dim(x))) stop('"x" must be a matrix-like object.')
  if((ncol(x) < 3L) || (ncol(x) > 4L)) stop('"x" must have 3 or 4 columns.')
  if(!isTRUE(all.equal(rowSums(x), rep(1.0, nrow(x)), check.attributes = FALSE))) stop('all values in "x" must be in [0, 1].')
  if(any((x < 0) | (x > 1))) stop('all values in "x" must be in [0, 1].')

  # transformations
  if(ncol(x) == 3L){
    return(toTernary(x))
  } else if(ncol(x) == 4L){
    return(toQuaternary(x))
  } else {
    stop("unexpected error")
  }
}
