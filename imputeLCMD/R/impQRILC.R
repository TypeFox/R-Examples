# -------------------------------------------------------------------------------------
# this function performs missing values imputation based quantile regression
# -------------------------------------------------------------------------------------

# arguments ___________________________________________________________________________
#           : dataSet.mvs - matrix of protein abundances with MVs 
#                                 (either peptides or proteins)
#           : tune.sigma - coefficient that controls the sd of the MNAR distribution
#                        - tune.sigma = 1 if the complete data distribution is supposed 
#                          to be gaussian
#                        - 0 < tune.sigma < 1 if the complete data distribution is  
#                          supposed to be left-censored

# output ______________________________________________________________________________
#           : results     - a list  containing: a matrix with the 
#                           complete abundances, a list with the 
#                           estimated parameters of the complete
#                           data distribution

impute.QRILC = function(dataSet.mvs,tune.sigma = 1){
  
  # get the dimension of the data .....................................................
  nFeatures = dim(dataSet.mvs)[1]
  nSamples = dim(dataSet.mvs)[2]
  
  # - initialize the matrix of complete abundances ......................................
  dataSet.imputed = dataSet.mvs
  
  # - initialize QR.obj which contains estimates of the distribution parameters
  #   for all samples
  QR.obj = list()
  
  for (i in 1:nSamples){
    
    curr.sample = dataSet.mvs[,i]
    
    # - calculate the percentage of missing values ........................................
    pNAs = length(which(is.na(curr.sample)))/length(curr.sample)
    
    # - estimate the mean and standard deviation of the original
    #   distribution using quantile regression
    
    upper.q = 0.99
        
    q.normal = qnorm(seq((pNAs+0.001),(upper.q+0.001),(upper.q-pNAs)/(upper.q*100)), 
                     mean = 0, sd = 1)
    
    q.curr.sample = quantile(curr.sample,
                             probs = seq(0.001,(upper.q+0.001),0.01),
                             na.rm = T) 

    temp.QR = lm(q.curr.sample ~ q.normal)
    QR.obj[[i]] = temp.QR
    
    mean.CDD = temp.QR$coefficients[1]
    sd.CDD = as.numeric(temp.QR$coefficients[2])
    
    # generate data from multivariate distributions with MLE parameters ....................
    data.to.imp = rtmvnorm(n=nFeatures, 
                           mean = mean.CDD, 
                           sigma = sd.CDD*tune.sigma,
                           upper = qnorm((pNAs+0.001),
                                         mean = mean.CDD,
                                         sd = sd.CDD),
                           algorithm=c("gibbs"))
    
    curr.sample.imputed = curr.sample
    curr.sample.imputed[which(is.na(curr.sample))] = data.to.imp[which(is.na(curr.sample))] 
    
    dataSet.imputed[,i] = curr.sample.imputed
    
  }
  
  results = list(dataSet.imputed,QR.obj)
  return(results)
  
}