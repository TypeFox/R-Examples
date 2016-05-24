oclKronVectMult1col <- function(kernel, a, b, y, ...){
  if ((ncol(a) * ncol(b)) != length(y)) stop("Error")
  matrix(oclRun(kernel, nrow(a) * nrow(b), as.double(as.vector(a)),
                as.double(as.vector(b)), as.double(y), ncol(a),
                nrow(a), ncol(b), nrow(b)),
         nrow = nrow(a) * nrow(b), byrow = TRUE)
}

