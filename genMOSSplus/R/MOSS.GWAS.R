MOSS.GWAS <-
function (alpha = 1, c = 0.1, cPrime = 0.0001, q = 0.1, replicates = 5, maxVars = 3, data, dimens, k = 2) {
  
  if (!is.data.frame(data)) {
    stop ("class(data) != 'data.frame'")
  }
  else if (length(dimens) != dim(data)[2]) {
    stop ("length(dimens) != dim(data)[2]")
  } 

  toKeep <- complete.cases(data) 
  data <- data[toKeep,,drop = F]
  dimens <- dimens[toKeep]

  checkArgsMOSS.GWAS (alpha, c, cPrime, q, replicates, data, maxVars, dimens, k)

  if (dim(data)[2] < maxVars)
    maxVars <- dim(data)[2]
  
  return (MOSS.GWAS.main(alpha, c, cPrime, q, replicates, maxVars, data, dimens, k))

}
