# -------------------------------------------------------------------------------------
# this function performs missing values imputation based on SVD algorithm
# -------------------------------------------------------------------------------------

# arguments ___________________________________________________________________________
#           : dataSet.mvs       - expression matrix with MVs (either peptides or proteins)
#           : K                 - the number of PCs

# output ______________________________________________________________________________
#           : dataSet.imputed   - expression matrix with MVs imputed

impute.wrapper.SVD = function(dataSet.mvs,K){
  
  resultSVD = pca(dataSet.mvs, method="svdImpute", nPcs = K)
  dataSet.imputed = resultSVD@completeObs
  
  return(dataSet.imputed)
  
}