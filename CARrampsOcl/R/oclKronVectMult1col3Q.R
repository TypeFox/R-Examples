oclKronVectMult1col3Q<- function(kernel, a, b, c, y, ...){
  if ((ncol(a) * ncol(b) * ncol(c)) != length(y)) stop("Error msg")
  matrix(oclRun(kernel, nrow(a) * nrow(b)*nrow(c), as.double(as.vector(a)),
                as.double(as.vector(b)), as.double(as.vector(c)),
                as.double(y), ncol(a),
                nrow(a), ncol(b), nrow(b), ncol(c), nrow(c)),
         nrow = nrow(a) * nrow(b)*nrow(c), byrow = TRUE)
}

