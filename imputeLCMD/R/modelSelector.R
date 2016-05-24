# .....................................................................................
# - this function determines row in the data matrix affected by a MNAR missingness 
#   mechanism
# - it is based on the assumption that the distributions of the mean values  
#   of proteins follows a normal distribution
# - the method makes use of a decision function defined as a tradeoff between the empirical 
#   CDF of the proteins' means and the theoretical CDF assuming that no MVs are present

# arguments ___________________________________________________________________________
#   : dataSet.mvs       - expression matrix containing abundances with MVs 
#                         (either peptides or proteins)

# output ______________________________________________________________________________
#   : results           - flags vector; "1" denotes rows containing random missing
#                         values; "0" denotes rows containing left-censored missing values

model.Selector = function(dataSet.mvs){
  
  # ___________________________________________________________________________________  
  # get the dimension of the data
  # -----------------------------------------------------------------------------------  
  nFeatures = dim(dataSet.mvs)[1]
  nSamples = dim(dataSet.mvs)[2] 
  
  # ___________________________________________________________________________________  
  # initialize the model selector variable
  # -----------------------------------------------------------------------------------  
  model.selector = rep(0,nFeatures)  
  
  # ___________________________________________________________________________________
  # compute the mean values for each protein
  # -----------------------------------------------------------------------------------
  mean.vect = rowMeans(dataSet.mvs, na.rm = T)
  
  # ___________________________________________________________________________________
  # estimate the percentage of missing values in the mean vector
  # -----------------------------------------------------------------------------------
  pNAs = length(which(is.na(mean.vect)))/length(mean.vect)
  
  # ___________________________________________________________________________________
  # - estimate the mean and standard deviation of the original
  #   distribution of the means using quantile regression
  # -----------------------------------------------------------------------------------
  
  upper.q = 0.99
  
  q.normal = qnorm(seq((pNAs+0.001),(upper.q+0.001),(upper.q-pNAs)/(upper.q*100)), 
                   mean = 0, sd = 1)
  
  q.mean.vect = quantile(mean.vect,
                           probs = seq(0.001,(upper.q+0.001),0.01),
                           na.rm = T) 
  
  temp.QR = lm(q.mean.vect ~ q.normal)
  QR.obj = temp.QR
  
  mean.CDD = temp.QR$coefficients[1]
  sd.CDD = as.numeric(temp.QR$coefficients[2])
  
  # _______________________________________________________________________________________
  # estimate the the censoring threshold
  # ---------------------------------------------------------------------------------------
  
  nPoints = 512
  
  # - create the support of the complete data distribution to evaluate the CDFs
  min.support = mean.CDD - 4*sd.CDD
  max.support = mean.CDD + 4*sd.CDD
  
  # - empirical estimation of the cdf of the observed data 
  
  mean.vect.sorted = sort(mean.vect)
  Fn <- ecdf(mean.vect.sorted)
  
  # - discretize the support of the CDFs
  
  support = c(seq(min.support,
                  min(mean.vect,na.rm = T),
                  length = nPoints),
              mean.vect.sorted,
              seq(max(mean.vect,na.rm = T),
                  max.support,
                  length = nPoints))
  
  # - evaluate the empirical cdf in the discrete points given by support
  
  cdf.OD = Fn(support)
  
  # - evaluate the theoretical cdf with MLE parameters estimated above 
  #   in the discrete points given by support
  
  cdf.CD = pnorm(support, 
                 mean = mean.CDD, 
                 sd = sd.CDD)
  
  # - compute the difference between the CDFs: evaluated in a particular point, 
  #   it gives the proportion of true discoveries (an estimation of the proportion
  #   of the true samples smaller than the evaluation point, given the data)
  
  diff.cdf = (cdf.CD - cdf.OD)
  
  # - the function to be optimize in order to obtain the optimum censoring threshold
  #   at its optimum, the function should maximize diff.cdf and should minimize 
  #   in the same time the cdf.OD
  
  # obj.fnc = diff.cdf-cdf.OD
  # obj.fnc = (diff.cdf)/(cdf.OD-1)
  # obj.fnc[is.infinite(obj.fnc)] = 0 
  # obj.fnc = 10*obj.fnc
  obj.fnc = (diff.cdf+1)/(cdf.OD+1) - 1
  obj.fnc = obj.fnc*10
  
  if (max(obj.fnc) > 0){
    censoring.thr = support[which(obj.fnc == max(obj.fnc))]
  }
  else{
    censoring.thr = min.support-1
  }
  
  # - set the model selector variable to "1" for the proteins with a mean value 
  #   higher than the KS statistic
  
  model.selector[which(mean.vect > censoring.thr)] = 1    
  
  result = list(model.selector,censoring.thr)
  
  return(result)
  
}