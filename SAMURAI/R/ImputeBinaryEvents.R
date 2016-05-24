ImputeBinaryEvents <-
function(table, effect.sizes.list, seed=NA, sims=1){
  # Imputes the number of events in the intervention arms of unpublished studies,
  # depending on the outlook(s) of those studies. 
  # 
  # Args: 
  #   table:  A data set with the following information:
  #             published studies, both arms: event counts, sample sizes
  #             unpublished studies:
  #               outlooks 
  #               control arm:      sample sizes, event counts
  #               intervention arm: sample sizes
  #   effect.sizes.list: A list of effect sizes to be assigned depending on outlook.  
  #   seed: random number seed
  #   sims: number of simulations (over which to average out the summary effect size)
  #
  # Returns:  The same data set with imputed event counts in intervention arms of unpublished studies.

  set.seed(NULL)  ## reset seed - redundant if done within a function call?
  if(!is.na(seed)){
    set.seed(seed)  
  }  

  # table <- Hpylori; sims <- 1       # testing

  ## Impute ctrl.events = proportion of events across published studies
  pub <- ExtractPublishedStudies(table)
  ctrl.propn <- sum(pub$ctrl.events) / sum(pub$ctrl.n)

  table <- ImputeControlGroupEvents(table)
  table <- CompleteOutlooksFactor(table)

  outlooks <- c("very positive", "positive", "no effect", "negative", "very negative", 
                "very positive CL", "positive CL", "current effect", "negative CL", "very negative CL")
  varnames <- c("vpos", "pos", "noef", "neg", "vneg", "vposcl", "poscl", "curr", "negcl", "vnegcl")

  # The number of different outlooks in the list above. 
  n <- length(outlooks)
   
  ## go through each of the outlooks in the list 
  for(i in 1:n){  
   
    # i <- 1  # testing

    # does any unpub in data set have that particular outlook?
    if(outlooks[i] %in% table$outlook){  
      
      # how many unpub studies have that outlook? 
      k <- length(which(table$outlook == outlooks[i]))  

      # vector of sample sizes of intervention arms
      expt.n <- table[which(table$outlook == outlooks[i]), ]$expt.n  

      # retrieve effect size for the outlook 
      # and multiply by the proportion of events across control arms
      # to impute the rate of events in the intervention arms
      expt.propn <- with(effect.sizes.list, get(varnames[i])) * ctrl.propn  
      
      # generate random numbers based on binomial distribution 
      # variable 'sims' determines how many simulations to be averaged over 
      table[which(table$outlook == outlooks[i]), ]$expt.events <- rbinom(n=k, size=expt.n*sims, prob=expt.propn) / sims  ## vector
    }  
  }  

  # Round to nearest integer
  table$ctrl.events <- round(table$ctrl.events)
  table$expt.events <- round(table$expt.events)
  # Avoid having # events exceed arm sample size # 
  table$ctrl.events <- pmin(table$ctrl.events, table$ctrl.n)
  table$expt.events <- pmin(table$expt.events, table$expt.n)
  # Avoid having # events negative
  table$ctrl.events <- pmax(table$ctrl.events, 0)
  table$expt.events <- pmax(table$expt.events, 0)
  
  return(table)
}
