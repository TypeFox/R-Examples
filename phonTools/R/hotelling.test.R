# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


hotelling.test <-
function (matrix1, matrix2 = NULL){

  if (is.null(ncol (matrix1))) return (cat ('Error: Univariate variable provided. Use a t-test.\n\n'))
  
  if (is.null(matrix2)){
    samples = 1
    mus1 = colMeans (matrix1); n1 = nrow (matrix1); p = ncol (matrix1)
    df1 = p; df2 = n1 - p; covar = solve (var (matrix1)/n1)
    f.value = t(mus1) %*% covar %*% (mus1) * ((n1 - p) / (p *(n1 -1)))
  }
  if (!is.null(matrix2)){
    samples = 2
    if (is.null(ncol (matrix2))) return (cat ('Error: Univariate variable provided. Use a t-test.\n\n'))
    if (ncol(matrix1) != ncol (matrix2)) return (cat ('Error: Variable dimensions do not match.\n\n'))

    mus1 = colMeans (matrix1); n1 = nrow (matrix1); p = ncol (matrix1); df1 = p;
    mus2 = colMeans (matrix2)
    n2 = nrow (matrix2)
    df2 = n1 + n2 - p - 1

    covar =  solve (((var(matrix1)*(n1-1) + var(matrix2)*(n2-1)) / (n1+n2-2)))
	
    f.value = t(mus1 - mus2) %*% covar %*% (mus1 - mus2) * ((n1*n2) / (n1+n2)) 
    f.value = f.value * ((n1 + n2 - p - 1) / (p*(n1 + n2 - 2)))  
  }
  print (paste(f.value,df1,df2))
  p.value = 1 - pf (f.value, df1, df2)

  output = list (f.value = f.value, df1 = df1, df2 = df2, p.value = p.value, samples = samples)
  class(output) = 'hotelling.test'
  output
}
