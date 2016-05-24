# -------------------------------------------------------------------------------------
# this function performs missing values imputation based on KNN algorithm
# -------------------------------------------------------------------------------------

# arguments ___________________________________________________________________________
#           : dataSet.mvs       - expression matrix with MVs (either peptides or proteins)
#           : K                 - the number of neighbors

# output ______________________________________________________________________________
#           : dataSet.imputed      - dataset containing complete abundances

impute.wrapper.KNN = function(dataSet.mvs,K){
  
  resultKNN = impute.knn(dataSet.mvs ,k = K, rowmax = 0.99, colmax = 0.99, maxp = 1500, rng.seed = sample(1:1000,1))
  dataSet.imputed = resultKNN[[1]]
  
  return(dataSet.imputed)
  
}