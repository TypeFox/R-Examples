HMM_simulation <-
function(size, m, delta = rep(1 / m, times = m), gamma = 0.8 * diag(m) + rep(0.2 / m, times = m), distribution_class, distribution_theta, obs_range=c(NA,NA), obs_round=FALSE, obs_non_neg = FALSE, plotting = 0)
{

######################################################################################################################################################################################################## Setting default-values for delta and gamma (if not available) ########################################################################################
################################################################################################################################################################################	
  if (!exists('delta')) 
  {
    delta <- rep(1 / m, times = m)
  }
  
  if (!exists('gamma')) 
  {
    gamma <- 0.8 * diag(m) + rep(0.2 / m, times = m)
  }
  
######################################################################################################################################################################################################## Generating of a random hidden Markov chain of the HMM(delta, gamma, distribution_theta, distribution_class) ##########################################
################################################################################################################################################################################
  
  markov_chain <- sample(x = seq(1, m, by = 1), 1, prob = delta)
  for (i in 2:size) 
  {
    last_state <- markov_chain[i - 1]
    markov_chain <- c(markov_chain, sample(x = seq(1, m, by=1), 1, prob = gamma[last_state, ]))
  }

######################################################################################################################################################################################################## Generating of random (hidden Markov chain)-dependend observations of the HMM(delta, gamma, distribution_theta, distribution_class) ###################
################################################################################################################################################################################  
  
  observations <- rep(NA, times = size)
   
  ######################### For a "pois"-HMM #################################################################################################################################
  if (distribution_class == "pois") 
  {	
    obs_dist_means <- distribution_theta$lambda
    
    means_along_markov_chain <- NULL
    for (i in 1:size) 
    {
      means_along_markov_chain <- c(means_along_markov_chain, distribution_theta$lambda[markov_chain[i]])	
    }
    
    for (i in 1:size)
    { 
      observations[i] <- rpois(n = 1, lambda = distribution_theta$lambda[markov_chain[i]])
      
      if (any(!is.na(obs_range))) 
      {
        if (!is.na(obs_range[1]) & !is.na(obs_range[2]))
        {
          while(observations[i] < obs_range[1] | observations[i] > obs_range[2])
          {
            observations[i] <- rpois(n = 1, lambda = distribution_theta$lambda[markov_chain[i]])
          }
        }
        if (!is.na(obs_range[1]) & is.na(obs_range[2]))
        {
          while(observations[i] < obs_range[1])
          {
            observations[i] <- rpois(n = 1, lambda = distribution_theta$lambda[markov_chain[i]])				 		
          }	
        }
        if (is.na(obs_range[1]) & is.na(obs_range[2]))
        {
          while(observations[i] > obs_range[2])
          {
            observations[i] <- rpois(n = 1, lambda = distribution_theta$lambda[markov_chain[i]])				 		
          }	
        }
      }
      
    }
  }
  
  
  ######################### For a "norm"-HMM #################################################################################################################################  
  if (distribution_class == "norm")
  {
    obs_dist_means <- distribution_theta$mean
    
    means_along_markov_chain <- NULL
    for(i in 1:size)
    {
      means_along_markov_chain <- c(means_along_markov_chain, distribution_theta$mean[markov_chain[i]])	
    }
    
    
    for (i in 1:size) 
    { 
      observations[i] <- rnorm(n = 1, mean = distribution_theta$mean[markov_chain[i]], sd = distribution_theta$sd[markov_chain[i]])
      
      if (any(!is.na(obs_range))) 
      {
        if (!is.na(obs_range[1]) & !is.na(obs_range[2])) 
        {
          while(observations[i] < obs_range[1] | observations[i] > obs_range[2])
          {
            observations[i] <- rnorm(n = 1, mean = distribution_theta$mean[markov_chain[i]], sd = distribution_theta$sd[markov_chain[i]])
          }
        }
        if (!is.na(obs_range[1]) & is.na(obs_range[2])) 
        {
          while(observations[i] < obs_range[1])
          {
            observations[i] <- rnorm(n = 1, mean = distribution_theta$mean[markov_chain[i]], sd = distribution_theta$sd[markov_chain[i]])
          }	
        }
        if (is.na(obs_range[1]) & is.na(obs_range[2])) 
        {
          while(observations[i] > obs_range[2])
          {
            observations[i] <- rnorm(n = 1, mean = distribution_theta$mean[markov_chain[i]], sd = distribution_theta$sd[markov_chain[i]])
          }	
        }
      }
      
      if (obs_non_neg == TRUE) 
      {
        if (observations[i] < 0)
        {
          observations[i] <- 0
        }	
      }
      
    }			
  }
  
  
  ######################### For a "genpois"-HMM ##############################################################################################################################  
  if (distribution_class == "genpois")
  {
    obs_dist_means <- distribution_theta$lambda1 / (1 - distribution_theta$lambda2)
    
    means_along_markov_chain <- NULL
    for (i in 1:size)
    {
      means_along_markov_chain <- c(means_along_markov_chain, (distribution_theta$lambda1[markov_chain[i]]) / (1 - distribution_theta$lambda2[markov_chain[i]]))	
    }
    
    
    for (i in 1:size)
    { 
      observations[i] <- rgenpois(n = 1, lambda1 = distribution_theta$lambda1[markov_chain[i]], lambda2 = distribution_theta$lambda2[markov_chain[i]])
      
      if (any(!is.na(obs_range)))
      {
        if (!is.na(obs_range[1]) & !is.na(obs_range[2]))
        {
          while(observations[i] < obs_range[1] | observations[i] > obs_range[2])
          {
            observations[i] <- rgenpois(n = 1, lambda1 = distribution_theta$lambda1[markov_chain[i]], lambda2 = distribution_theta$lambda2[markov_chain[i]])
          }
        }
        if (!is.na(obs_range[1]) & is.na(obs_range[2]))
        {
          while(observations[i] < obs_range[1])
          {
            observations[i] <- rgenpois(n = 1, lambda1 = distribution_theta$lambda1[markov_chain[i]], lambda2 = distribution_theta$lambda2[markov_chain[i]])				 		}	
        }
        if(is.na(obs_range[1]) & is.na(obs_range[2]))
        {
          while(observations[i] > obs_range[2])
          {
            observations[i] <- rgenpois(n = 1, lambda1 = distribution_theta$lambda1[markov_chain[i]], lambda2 = distribution_theta$lambda2[markov_chain[i]])
          }	
        }
      }				
      
    }
  }
  
  ######################### For a "geom"-HMM #################################################################################################################################
  if (distribution_class == "geom")
  {	
    obs_dist_means <- distribution_theta$prob
    
    means_along_markov_chain <- NULL
    for (i in 1:size)
    {
      means_along_markov_chain <- c(means_along_markov_chain, distribution_theta$prob[markov_chain[i]])	
    }
    
    for (i in 1:size)
    { 
      observations[i] <- rgeom(n = 1, prob = distribution_theta$prob[markov_chain[i]])
      
      if (any(!is.na(obs_range))) 
      {
        if (!is.na(obs_range[1]) & !is.na(obs_range[2])) 
        {
          while(observations[i] < obs_range[1] | observations[i] > obs_range[2]) 
          {
            observations[i] <- rgeom(n = 1,prob = distribution_theta$prob[markov_chain[i]]) 
          }
        }
        if (!is.na(obs_range[1]) & is.na(obs_range[2])) 
        {
          while(observations[i] < obs_range[1])
          {
            observations[i] <- rgeom(n = 1,prob = distribution_theta$prob[markov_chain[i]])
          }	
        }
        if (is.na(obs_range[1]) & is.na(obs_range[2]))
        {
          while(observations[i] > obs_range[2])
          {
            observations[i] <- rgeom(n = 1, prob = distribution_theta$prob[markov_chain[i]])
          }	
        }
      }
    }
  }
  
  
  
######################################################################################################################################################################################################## Plotting of the outcome ##############################################################################################################################
################################################################################################################################################################################  

  if (!is.na(plotting)) 
  {	
    
    if (plotting == 0)
    {   
      par(mfrow=c(2,2))	
      
      plot(markov_chain, xlab ='t', main ='simulated (hidden) Markov chain', col = "green", type = 'o', ylab = "states")
      
      plot(means_along_markov_chain, xlab = 't', main = 'means along Markov chain', col = "green", type = 'o', ylab = 'observation')
      
      plot(observations, xlab = 't', main = 'observations along Markov chain')
      abline(h = obs_dist_means, col = "grey50", lty = "dashed")
      lines(means_along_markov_chain, xlab = 'time', main = 'simulation data', type = "l", col = "green")
      
      plot(observations, xlab = 't', main = 'simulated observations')
      
      par(mfrow=c(1,1))
      
      plot(markov_chain, xlab = 't', main = 'simulated (hidden) Markov chain', col = "green", type = 'o', ylab = "states")
      
      
      plot(means_along_markov_chain, xlab = 't', main = 'means along Markov chain', col = "green", type = 'o', ylab = 'observation')
      
      plot(observations, xlab = 't', main = 'observations along Markov chain')
      abline(h = obs_dist_means, col = "grey50", lty = "dashed")
      lines(means_along_markov_chain, xlab = 'time', main = 'simulation data', type = "l", col = "green")
      
      plot(observations, xlab = 't', main = 'simulated observations')
      
      par(mfrow=c(1,1))  
      
    }
    if (plotting == 1) 
    {   
      par(mfrow=c(2,2))	
      
      plot(markov_chain, xlab = 't', main = 'simulated (hidden) Markov chain', col = "green", type = 'o', ylab = "states")
      
      
      plot(means_along_markov_chain, xlab = 't', main = 'means along Markov chain', col = "green", type = 'o', ylab = 'observation')
      
      plot(observations, xlab = 't', main = 'observations along Markov chain')
      abline(h = obs_dist_means, col = "grey50", lty = "dashed")
      lines(means_along_markov_chain, xlab = 'time', main = 'simulation data', type = "l", col = "green")
      
      plot(observations, xlab = 't', main = 'simulated observations')
      
      par(mfrow=c(1,1))
    }
    
    if (plotting == 2) 
    {   
      plot(markov_chain, xlab = 't', main = 'simulated (hidden) Markov chain', col = "green", type = 'o', ylab = "states")
    }
    if (plotting == 3)
    {   
      plot(means_along_markov_chain, xlab = 't', main = 'means along Markov chain', col = "green", type = 'o', ylab = 'observation')
    }
    if (plotting == 4) 
    {   
      plot(observations, xlab = 't', main = 'observations along Markov chain')
      abline(h = obs_dist_means, col = "grey50", lty = "dashed")
      lines(means_along_markov_chain, xlab = 'time', main = 'simulation data', type = "l", col = "green")
    }
    if (plotting == 5)
    {   
      plot(observations, xlab = 't', main = 'simulated observations')
    }
  }
  
  if(obs_round == TRUE)
  {
    observations <- round(observations)	
  }	

############################################################################################################################################################################################################################### Return results ################################################################################################################
################################################################################################################################################################################
  
  return(list(size = size, 
  			      m = m, 
  			      delta = delta, 
  			      gamma = gamma, 
  			      distribution_class = distribution_class, 
  			      distribution_theta = distribution_theta, 
  			      markov_chain = markov_chain, 
  			      means_along_markov_chain = means_along_markov_chain, 
  			      observations = observations))
}             