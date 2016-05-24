# .....................................................................................
# this function performs missing values imputation by 0

# arguments ___________________________________________________________________________
#           : dataSet.mvs       - MSnSet containing abundances with MVs (either peptides or proteins)

# output ______________________________________________________________________________
#           : dataSet.imputed      - dataset containing complete abundances

impute.ZERO = function(dataSet.mvs){
  
  dataSet.imputed = dataSet.mvs
  dataSet.imputed[which(is.na(dataSet.mvs))] = 0
  
  return(dataSet.imputed)
  
}