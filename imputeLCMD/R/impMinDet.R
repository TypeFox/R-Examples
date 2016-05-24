# ---------------------------------------------------------------------------------
# this function performs missing values imputation by the minimum value observed

# input ___________________________________________________________________________
#         : dataSet.mvs       - expression matrix containing abundances with 
#                               MVs (either peptides or proteins);
#         : q                 - the q quantile used to estimate the minimum 
#                               value observed for each sample; 

# output __________________________________________________________________________
#         : dataSet.imputed      - dataset containing complete abundances

impute.MinDet = function(dataSet.mvs,q = 0.01){
  
  nSamples = dim(dataSet.mvs)[2]
  dataSet.imputed = dataSet.mvs
  
  lowQuantile.samples = apply(dataSet.imputed,2,quantile,prob = q,na.rm = T)
  
  for (i in 1:(nSamples)){
    dataSet.imputed[which(is.na(dataSet.mvs[,i])),i] = lowQuantile.samples[i]
  }
    
  return(dataSet.imputed)
  
}