# Create an orthogonal basis of K functions in [0, 1], with nGrid points.
# Output: a K by nGrid matrix, each column containing an basis function.
CreateBasis <- function(K, pts=seq(0, 1, length.out=50), type='sin') {

  nGrid <- length(pts)
  possibleTypes <- c('cos', 'sin', 'fourier', 'unknown')
  type <- possibleTypes[pmatch(type, possibleTypes, nomatch=length(possibleTypes))]
  
  stopifnot(is.numeric(K) && length(K) == 1 && K > 0)
  
  if (type == 'cos') {
    sapply(seq_len(K), function(k) 
      if (k == 1) {
        rep(1, nGrid)
      } else {
        sqrt(2) * cos((k - 1) * pi * pts)
      }
    )
  } else if (type == 'sin') {
    sapply(seq_len(K), function(k) sqrt(2) * sin(k * pi * pts))
  } else if (type == 'fourier') {
    sapply(seq_len(K), function(k) 
      if (k == 1) {
        rep(1, nGrid)
      } else if (k %% 2 == 0) {
        sqrt(2) * sin(k * pi * pts)
      } else {
        sqrt(2) * cos((k - 1) * pi * pts)
      }
    )
  } else if (type == 'unknown') {
    stop('unknown basis type')
  }
}