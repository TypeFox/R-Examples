################################################################################
################################################################################
## author Till Junge <till.junge@altermail.ch>                                ##
##                                                                            ##
## Copyright (c) UNIL (Universite de Lausanne)                                ##
## NCCR - LIVES (National Centre of Competence in Research "LIVES -           ##
## Overcoming vulnerability: life course perspectives",                       ##
## <http://www.lives-nccr.ch/>)                                               ##
##                                                                            ##
## spacom is free software: you can redistribute it and/or modify it under    ##
## the terms of the GNU General Public License as published by the Free       ##
## Software Foundation, either version 2 of the License or any later version. ##
##                                                                            ##
## spacom is distributed in the hope that it will be useful, but WITHOUT ANY  ##
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  ##
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more      ##
## details, see <http://www.gnu.org/licenses/>.                               ##
################################################################################
################################################################################

makeWeightsObject <- function(distance.matrix,
                              kernel,
                              moran) {
  obj <- new("weightsObject")
  if (is(distance.matrix, "data.frame")) {
    distance.matrix <- as.matrix(distance.matrix)
  }
  if (!is(distance.matrix, "matrix") && !is(distance.matrix, "Matrix")) {
    stop("The distance matrix has to be of class 'matrix' or 'Matrix'. You ",
         "specified an object of class '", class(distance.matrix), "'.")
  }

  if (!nrow(distance.matrix) == ncol(distance.matrix)) {
    stop("The distance matrix has to be square, you specified a matrix of ",
         "size ", nrow(distance.matrix), "x", ncol(distance.matrix))
  }

  obj@distance.matrix <- Matrix(distance.matrix)

  obj@kernel <- checkKernel(kernel)

  ##
  tryCatch(moran <- as.logical(moran),
           error=function(er) {
             stop("The variable you specified for the argument 'moran' could ",
                  "not be coerced into a logical value. Please specify TRUE ",
                  "or False")})
  obj@moran <- moran

  return(obj)
}

performWeights <- function(obj, bandwidth) {
  mat <- obj@kernel(obj@distance.matrix, bandwidth)
  if (obj@moran) {
    diag(mat) <- 0
  }

  return(mat)
}

WeightMatrix <- function(distance.matrix, bandwidth, kernel=NULL, moran=FALSE) {
  obj <- makeWeightsObject(distance.matrix, kernel, moran)
  return(performWeights(obj, bandwidth))
}
