MOSS_GWAS <-
function (alpha = 1, c = 0.1, cPrime = 0.0001, q = 0.1, replicates = 5, maxVars = 3, data, dimens, confVars = NULL, k = NULL) {
  
  if (!is.data.frame(data)) {
    stop ("class(data) != 'data.frame'")
  }
  else if (is.null(colnames(data))) {
    stop ("Please add column names to data")
  }
  
  toKeep <- complete.cases(data) 
  data <- data[toKeep,,drop = F]
  
  check_args_MOSS_GWAS (alpha, c, cPrime, q, replicates, maxVars, data, dimens, confVars, k)

  return (MOSS_GWAS_main(alpha, c, cPrime, q, replicates, maxVars, data, dimens, confVars, k))

}
