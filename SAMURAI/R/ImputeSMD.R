ImputeSMD <-
function(table, effect.sizes.list, seed=NA, noise=0.01){
  # Imputes the standardized mean difference of unpublished studies,
  # depending on the outlook(s) of those studies. 
  # 
  # Args: 
  #   table:  A data set with the following information:
  #             published studies: means, sd's, sample sizes
  #             unpublished studies: outlooks, sample sizes
  #   effect.sizes.list: A list of effect sizes to be assigned depending on outlook.  
  #   seed: random number seed
  #   noise: user added Gaussian random noise; standard deviation from the assigned effect size 
  #
  # Returns:  The same data set with imputed effect sizes for unpublished studies.

  set.seed(NULL)  ## reset seed - redundant if done within a function call?
  if(!is.na(seed)){
    set.seed(seed)  
  }  

  table <- CompleteOutlooksFactor(table)
  
  outlooks <- c("very positive", "positive", "no effect", "negative", "very negative", 
                "very positive CL", "positive CL", "current effect", "negative CL", "very negative CL")
  varnames <- c("vpos", "pos", "noef", "neg", "vneg", "vposcl", "poscl", "curr", "negcl", "vnegcl")
  n <- length(outlooks)
  
  for(i in 1:n){
    # i <- 4
    if(outlooks[i] %in% table$outlook){ # if there unpub assigned to that outlook
      k <- length(which(table$outlook == outlooks[i]))
      # retrieve smd to be assigned
      smd.assigned <- with(effect.sizes.list, get(varnames[i]))
      # add some random noise
      table[which(table$outlook==outlooks[i]), ]$yi <- rnorm(k, mean=smd.assigned, sd=noise)
    }  
  }  
  
  return(table)
}
