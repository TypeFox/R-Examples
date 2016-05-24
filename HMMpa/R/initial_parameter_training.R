initial_parameter_training <-
function(x, m, distribution_class, n = 100, discr_logL =  FALSE, discr_logL_eps = 0.5)
{	
  ################################################################################################################################################################################################################################# Setting fix values for all distributions ##################################################################################
##############################################################################################################################################################################
  
  par_var = 0.5
  gamma <- 0.8 * diag(m) + rep(0.2 / m, m)
  delta <- rep(1 / m, m)
  
  if (distribution_class == "pois" | distribution_class == "geom")
  {
    k <- 1
  }	
  if (distribution_class == "norm" | distribution_class == "genpois")
  {
    k <- 2
  }	
  
################################################################################################################################################################################################################################# Finding plausible starting values for "pois" ###############################################################################
################################################################################################################################################################################
  
  if (distribution_class == "pois") 
  {	
  	
    ############### Building array-object to save 2+n different sets of parameters #########################################################################################
    erg <- array(NA, dim = c(m, 1 + 1 + k, (n + 2)))    
    colnames(erg)<-c("logL", "E","lambda")
    
    
    
    ############### Bilding the first set of parameters via 'quantile'-based picking for the mean values E, and saving them into the first place of the array ##############
    E <- quantile(x = x, probs = seq(0, 1, 1 / (m)), na.rm = FALSE, names = FALSE, type = 7)
    E=E[1:m]
    E <- sort(E)
    distribution_theta <- list(lambda = E)
    fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
    if (inherits(fb,"try-error"))
    {  
      fb$logL <- -Inf
    }
    logL <- fb$logL
    erg[1,1,1] <- logL
    for (j in 1:m)
    {
      erg[j,2,1] <- E[j]
    }
    for (j in 1:m)
    {
      erg[j,3,1] <- distribution_theta$lambda[j]
    }
    
    
    
    ############## Bilding the second set of parameters via 'some kind of uniformly-distribution'-picked values as E, and saving them into the second place of the array ####
    E <- seq((min(x)+1), max(x), length = m)
    E <- sort(E)
    distribution_theta <- list(lambda = E)
    fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
    if (inherits(fb,"try-error"))
    {  
      fb$logL <- -Inf
    }
    logL <- fb$logL
    erg[1,1,2] <- logL
    for (j in 1:m)
    {
      erg[j,2,2] <- E[j]
    }
    for (j in 1:m)
    {
      erg[j,3,2] <- distribution_theta$lambda[j]
    }
    
    
    ############## Bilding the 3:n + 2 (loop) sets of parameters via 'randomly'-picked values as E, and saving them into the 3:n + 2 places of the array ########################
    if ( n > 0)
    {
      for (i in 3:(n + 2))
      {  
        xt <- x
        E <- sample(x = x, 1, replace = FALSE, prob = NULL)
        for (h in 1:(m - 1))
        {	
          xt <- xt[!xt == E[h]]
          E <- append(E, sample(x = xt, 1, replace = FALSE, prob = NULL))
        }
        E <- sort(E)
        distribution_theta <- list(lambda = E)
        fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta,discr_logL =  discr_logL, discr_logL_eps = 0.5), silent = FALSE)	
        if (inherits(fb,"try-error"))
        {  
          fb$logL <- -Inf
        }       
        logL <- fb$logL
        erg[1,1,i] <- logL		
        for (j in 1:m)
        {
          erg[j,2,i] <- E[j]
        }
        for (j in 1:m)
        {
          erg[j,3,i] <- distribution_theta$lambda[j]
        }
      }
    }	
    L_min <- which.max(erg[1,1,])
    return(list(m = m,k = k,logL = erg[1,1,L_min], E = erg[,2,L_min], distribution_theta = list(lambda = erg[,3,L_min]), delta = delta, gamma = gamma))	
  }
  
    
  
  
################################################################################################################################################################################################################################# Finding plausible starting values for "geom" ###############################################################################
################################################################################################################################################################################  
    
  if (distribution_class == "geom") 
  {	

    ############### Building array-object to save 2+n different sets of parameters #########################################################################################
    erg <- array(NA, dim = c(m, 1 + 1 + k, (n + 2)))    
    colnames(erg)<-c("logL", "E", "prob")



    ############### Bilding the first set of parameters via 'quantile'-based picking for the mean values E, and saving them into the first place of the array ##############
    E <- quantile(x = x, probs = seq(0, 1, 1 / (m)), na.rm = FALSE,names = F, type = 7)
    E=E[1:m]
    E <- sort(E)
    prob <- 1 / E
    distribution_theta <- list(prob = prob)
    fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
    if (inherits(fb,"try-error"))
    {  
      fb$logL <- -Inf
    }
    logL <- fb$logL
    erg[1,1,1] <- logL
    for (j in 1:m)
    {
      erg[j,2,1] <- E[j]
    }
    for (j in 1:m)
    {
      erg[j,3,1] <- distribution_theta$prob[j]
    }



    ############## Bilding the second set of parameters via 'some kind of uniformly-distribution'-picked values as E, and saving them into the second place of the array ####
    E <- seq((min(x) + 1), max(x), length = m)
    E <- sort(E)
    prob <- 1 / E
    distribution_theta <- list(prob = prob)
    fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
    if (inherits(fb,"try-error"))
    {  
      fb$logL <- -Inf
    }
    logL <- fb$logL
    erg[1,1,2] <- logL
    for (j in 1:m)
    {
      erg[j,2,2] <- E[j]
    }
    for (j in 1:m)
    {
      erg[j,3,2] <- distribution_theta$prob[j]
    }



    ############## Bilding the 3:n + 2 (loop) sets of parameters via 'randomly'-picked values as E, and saving them into the 3:n + 2 places of the array ########################
    if ( n > 0)
    {
      for (i in 3:(n + 2))
      {  
        xt <- x
        E <- sample(x = x, 1, replace = FALSE, prob = NULL)
        for (h in 1:(m-1))
        {	
          xt <- xt[!xt==E[h]]
          E <- append(E, sample(x = xt, 1,replace = FALSE, prob = NULL))
        }
        E <- sort(E)
        prob <- 1 / E
        distribution_theta <- list(prob = prob)
        fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta,discr_logL =  discr_logL,discr_logL_eps = 0.5), silent = FALSE)	
        if (inherits(fb,"try-error"))
        {  
          fb$logL <- -Inf
        }
        logL <- fb$logL
        erg[1,1,i] <- logL		
        for (j in 1:m)
        {
          erg[j,2,i] <- E[j]
        } 
        for (j in 1:m)
        {
          erg[j,3,i] <- distribution_theta$prob[j]
        }
      }
    }
    L_min <- which.max(erg[1,1,])
    return(list(m = m,k = k,logL = erg[1,1,L_min], E=erg[,2,L_min], distribution_theta = list(prob=erg[,3,L_min]), delta = delta, gamma = gamma))	
  }
    
  
  
  
  
################################################################################################################################################################################################################################# Finding plausible starting values for "norm" ###############################################################################
################################################################################################################################################################################    
    
  
  if (distribution_class == "norm") 
  {	
    
    
    ############### Building array-object to save 2+n different sets of parameters #########################################################################################
    erg <- array(NA, dim = c(m, 1 + 1 + k, (n + 2)))    
    colnames(erg)<-c("logL", "E", "mean", "sd")
   
   
   
    ############### Bilding the first set of parameters via 'quantile'-based picking for the mean values E, and saving them into the first place of the array ##############   
    E <- quantile(x = x, probs = seq(0, 1, 1/(m)), na.rm = FALSE,names = F, type = 7)
    E <- E[1:m]
    E <- sort(E)
    sd_obs_pi <- rep(sd(x) / m, times = m)
    distribution_theta <- list(mean = E, sd = sd_obs_pi)
    fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
    if (inherits(fb,"try-error"))
    {  
      fb$logL <- -Inf
    }
    logL <- fb$logL
    erg[1,1,1] <- logL
    for (j in 1:m)
    {
      erg[j,2,1] <- E[j]
    }
    for (j in 1:m)
    {
      erg[j,3,1] <- distribution_theta$mean[j]
    }
    for (j in 1:m)
    {
      erg[j,4,1] <- distribution_theta$sd[j]
    }
 
 
 
    ############## Bilding the second set of parameters via 'some kind of uniformly-distribution'-picked values as E, and saving them into the second place of the array #### 
    E <- seq((min(x) + 1), max(x), length = m)
    E <- sort(E)
    sd_obs_pi <- rep(sd(x) / m, times = m)
    distribution_theta <- list(mean = E, sd = sd_obs_pi)
    fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
    if (inherits(fb,"try-error"))
    {  
      fb$logL <- -Inf
    }
    logL <- fb$logL
    erg[1,1,2] <- logL
    for (j in 1:m)
    {
      erg[j,2,2] <- E[j]
    }
    for (j in 1:m)
    {
      erg[j,3,2] <- distribution_theta$mean[j]
    }
    for (j in 1:m)
    {
      erg[j,4,2] <- distribution_theta$sd[j]
    }




    ############## Bilding the 3:n + 2 (loop) sets of parameters via 'randomly'-picked values as E, and saving them into the 3:n + 2 places of the array ########################
    if ( n > 0)
    {
      for (i in 3:(n + 2))
      {   
        xt <- x
        E <- sample(x = x, 1, replace = FALSE, prob = NULL)
        for (h in 1:(m - 1))
        {	
          xt <- xt[!xt == E[h]]
          E <- append(E, sample(x = xt, 1,replace = FALSE, prob = NULL))
        }
        E <- sort(E)		
        sd_obs_pi <- rep(sd(x) / m, times = m)
        distribution_theta <- list(mean = E, sd = sd_obs_pi)
        fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
        if (inherits(fb,"try-error"))
        {  
          fb$logL <- -Inf
        }
        logL <- fb$logL
        erg[1,1,i] <- logL
        for (j in 1:m)
        {
          erg[j,2,i] <- E[j]
        }
        for (j in 1:m)
        {
          erg[j,3,i] <- distribution_theta$mean[j]
        }	
        for (j in 1:m)
        {
          erg[j,4,i] <- distribution_theta$sd[j]
        }	
      }
    }		
    L_min <- which.max(erg[1,1,])
    return(list(m = m,k = k,logL = erg[1,1,L_min], E=erg[,2,L_min], distribution_theta = list(mean = erg[,3,L_min], sd = erg[,4,L_min]), delta = delta, gamma = gamma))
  }
  

################################################################################################################################################################################################################################# Finding plausible starting values for "genpois" ############################################################################
################################################################################################################################################################################     
  
  if (distribution_class == "genpois") 
  {	

    ############### Building array-object to save 2+n different sets of parameters #########################################################################################
    erg <- array(NA, dim = c(m, 1 + 1 + k, (n + 2)))    
    colnames(erg)<-c("logL", "E","lambda1","lambda2")
 
 
    ############### Bilding the first set of parameters via 'quantile'-based picking for the mean values E, and saving them into the first place of the array ############## 
    E <- quantile(x = x, probs = seq(0, 1, 1/(m)), na.rm = FALSE,names = F, type = 7)
    E <- E[1:m]
    E <- sort(E)
    E[E == 0] <- 1
    lambda2 <- rep(par_var, times = m)
    lambda1 <- E * (1 - lambda2)
    distribution_theta <- list(lambda1 = lambda1, lambda2 = lambda2)
    fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
    if (inherits(fb,"try-error"))
    {  
      fb$logL <- -Inf
    }
    logL <- fb$logL
    erg[1,1,1] <- logL
    for (j in 1:m)
    {
      erg[j,2,1] <- E[j]
    }
    for (j in 1:m)
    {
      erg[j,3,1] <- distribution_theta$lambda1[j]
    }
    for (j in 1:m)
    {
      erg[j,4,1] <- distribution_theta$lambda2[j]
    }



    ############## Bilding the second set of parameters via 'some kind of uniformly-distribution'-picked values as E, and saving them into the second place of the array ####
    E <- seq((min(x) + 1), max(x), length = m)
    E <- sort(E)
    E[E == 0] <- 1
    lambda2 <- rep(par_var, times = m)
    lambda1 <- E * (1 - lambda2)
    distribution_theta <- list(lambda1 = lambda1, lambda2 = lambda2)
    fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
    if (inherits(fb,"try-error"))
    {  
      fb$logL <- -Inf
    }
    logL <- fb$logL
    erg[1,1,2] <- logL
    for (j in 1:m)
    {
      erg[j,2,2] <- E[j]
    }
    for (j in 1:m)
    {
      erg[j,3,2] <- distribution_theta$lambda1[j]
    }
    for (j in 1:m)
    {
      erg[j,4,2] <- distribution_theta$lambda2[j]
    }



    ############## Bilding the 3:n + 2 (loop) sets of parameters via 'randomly'-picked values as E, and saving them into the 3:n + 2 places of the array ########################
    if ( n > 0)
    {
      for (i in 3:(n + 2))
      {  
        xt <- x
        E=sample(x = x, 1,replace = FALSE, prob = NULL)
        for (h in 1:(m - 1))
        {	
          xt <- xt[!xt == E[h]]
          E <- append(E, sample(x = xt, 1, replace = FALSE, prob = NULL))
        }
        E <- sort(E)
        E[E == 0] <- 1
        lambda2 <- rep(par_var, times = m)
        lambda1 <- E * (1 - lambda2)
        distribution_theta <- list(lambda1 = lambda1, lambda2 = lambda2)
        fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
        if (inherits(fb,"try-error"))
        {  
          fb$logL <- -Inf
        }       
        logL <- fb$logL
        erg[1,1,i] <- logL
        for (j in 1:m)
        {
          erg[j,2,i] <- E[j]
        }
        for (j in 1:m)
        {
          erg[j,3,i] <- distribution_theta$lambda1[j]
        }	
        for (j in 1:m)
        {
          erg[j,4,i] <- distribution_theta$lambda2[j]
        }		
      }
    }
    L_min <- which.max(erg[1,1,])
    return(list(m = m, k = k, logL = erg[1,1,L_min], E=erg[,2,L_min], distribution_theta=list(lambda1=erg[,3,L_min], lambda2=erg[,4,L_min]), delta = delta, gamma = gamma))
  }
  
  
}
