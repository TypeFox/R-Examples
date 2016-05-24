# -------------------------------------------------------------------------------------
# this function generates missing data in a complete data matrix
# -------------------------------------------------------------------------------------

# arguments ___________________________________________________________________________
#           : original        - complete data matrix containing all measurements
#
#           : mean.THR,sd.THR   - parameters of the threshold distribution which 
#                               controls the MVs rate (mean.THR should be initially set 
#                               such that the result of the initial thresholding, 
#                               in terms of no. of NAs, equals the desired total 
#                               missing data rate) 
#                             - example: if one wants to generate 30% missing data
#                               mean.THR can be set as follows: 
#                               mean.THR = quantile(pepExprsData, probs = 0.3)
#                             - sd.THR is usually set to a small value (e.g. 0.1)
#
#           : MNAR.rate       - percentage of MVs which are missing not at random
#
# output ______________________________________________________________________________
#           : list            - contains the original complete data matrix, 
#                               the data matrix with missing data and the 
#                               percentage of missing data

insertMVs = function(original,mean.THR,sd.THR,MNAR.rate){
  
  originalNaNs = original
  
  nProt = nrow(original)
  nSamples = ncol(original)
  
  # -------------------------------------------------------------------------------
  # insert MNAR-MVs 
  
  # - define threshold matrix
  thr = matrix(rnorm(nSamples * nProt, mean.THR, sd.THR), nProt, nSamples)
      
  # - find location of MNAR values
  indices.MNAR = which(original<thr)
  
  # - number of values to be replaced with MNARs
  no.MNAR = round(MNAR.rate/100*length(indices.MNAR))
      
  # insert a percentage of MVs which equals (100-MNAR.rate(j))*length(indices_MNAR)
  temp = matrix(original,1,nSamples * nProt)
  temp[sample(indices.MNAR,no.MNAR)] = NA
      
  
  # -------------------------------------------------------------------------------
  # insert MCAR-MVs with a rate of MCAR_rate(j)
  
  # - find location of potential MCAR values
  indices.MCAR = which(!is.na(temp))
  
  # - number of values to be replaced with MCARs
  no.MCAR = floor((100-MNAR.rate)/100*length(indices.MNAR))
  print(no.MCAR+no.MNAR)
  
  # - insert a percentage of MVs which equals MCAR_rate(j)*length(indices_MNAR)
  temp[sample(indices.MCAR,no.MCAR)] = NA
      
  originalNaNs = matrix(temp,nProt,nSamples)
  
  # -------------------------------------------------------------------------------
  # adjust for variables in which all samples are 'NA'
  
  originalNaNs_adjusted = originalNaNs
  
  # ..............................................................
  noNaNs_Var = rowSums(is.na(originalNaNs))
  allNaNs_Vars = which(noNaNs_Var==nSamples)
  
  sampleIndexToReplace = sample(1:nSamples,length(allNaNs_Vars),replace = T)
  
  for (i in 0:length(sampleIndexToReplace))
  {
    originalNaNs_adjusted[allNaNs_Vars[i],sampleIndexToReplace[i]] = original[allNaNs_Vars[i],sampleIndexToReplace[i]]
  }
  # ..............................................................
  
  original.mvs = originalNaNs_adjusted
  
  # calculate the real percentage of MVs
  pNaNs = length(which(is.na(original.mvs)))/(nSamples*nProt)
  
  return(list(original,original.mvs,pNaNs))
  
}