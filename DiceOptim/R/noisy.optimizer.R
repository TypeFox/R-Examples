#source("noisy.EGO.R")
noisy.optimizer <- function(optim.crit, optim.param=NULL, model, n.ite, noise.var=NULL, funnoise, lower, upper, 
                      parinit=NULL, control=NULL,CovReEstimate=TRUE,NoiseReEstimate=FALSE, 
                      nugget.LB=1e-5, estim.model=NULL, type="UK")
{
  ############################################################################################################
  # A switch is made at every iteration to choose the corresponding infill criterion, then the model is updated 
  # similarly for all methods.
  #
  # All methods return a final noisy kriging model. In addition, the quantities "all.X", "all.y", "all.thetas"
  # (size=n.ite) are saved for analysis.
  #
  # This version allows sequential re-estimation of the noise variance. An additional model is created specifically
  # for parameter estimation (estim.model), with uniform nugget (repeated experiments are not aggregated).
  # "model" and "estim.model" are equivalent, but "model" has less observations ( hence goes faster).
  #
  # Additional (optional) arguments:
  #    - "NoiseReEstimate": TRUE/FALSE
  #    - "noise.var": initial guess for the unknown variance (used in optimization)
  #    - "obs.n.rep": number of repetitions per observation point. Required if "model" has heterogeneous variances
  #    - "estim.model": model with homogeneous nugget effect (no noise.var). Required only for restarting algorithms
  #
  ############################################################################################################
  
  cat("Starting noisy optimization with the following criterion and parameters \n")
  cat(optim.crit, "\n")
  if(!is.null(optim.param))  print(optim.param,quote=FALSE)
  cat("----------------- \n")
  optim.result <- list()
  
  #---------------------------------------------------------------------------------------------------------
  # Initialization for the unknown noise variance case
  #---------------------------------------------------------------------------------------------------------
  if (NoiseReEstimate)
  { if (!CovReEstimate)
  { cat("Noise variance cannot be re-estimated without re-estimating the covariance parameters \n")
    cat("covReEstimate switched to TRUE \n")
    CovReEstimate = TRUE
  }
    
    obs.n.rep <- pmax( 1, round(noise.var / as.numeric(model@noise.var)) ) 
    
    if (!is.null(estim.model))
      #-- Case 1.1: estim.model is provided (noise.var is set to the nugget of estim.model) -----------------
      #-- All the data is assumed to be in the appropriate format -------------------------------------------
    { noise.var <- estim.model@covariance@nugget
    } else
      
      #-- Case 1.2: estim.model needs to be built (noise.var is estimated) ----------------------------------
    { if ( max(obs.n.rep)==1 )
      { estim.model <- km(formula=model@trend.formula, covtype=model@covariance@name,
                         design=model@X, response=model@y,
                         lower=model@lower, upper=model@upper,
                          coef.cov=covparam2vect(model@covariance), coef.var=model@covariance@sd2,
                          optim.method=model@optim.method, 
                         nugget=noise.var, control=model@control)
        if (type=="SK"){ estim.model@trend.coef=model@trend.coef }
        estim.model@covariance@nugget.estim =TRUE
      } else 
      { cat("Unless estim.model is provided, noisy.EGO cannot estimate the noise variance when the initial 
            model has heterogeneous noise variance. NoiseReEstimate switched to FALSE. \n")
        NoiseReEstimate <- FALSE
      }
    }
  }

  ##########################################################################################################
  ### MAIN LOOP STARTS #####################################################################################
  for (i.time.steps in 1:n.ite)
  {
    add.obs <- TRUE

    # Update number of repetitions in the unknown noise variance case
    if (NoiseReEstimate) {  obs.n.rep <- pmax( 1, round(noise.var / as.numeric(model@noise.var)) )  }

    #-------------------------------------------------------------------------------------------------------
    #-- Find the new observation based on the REINTERPOLATION criterion ------------------------------------
    #-------------------------------------------------------------------------------------------------------
    if (optim.crit == "reinterpolation")
    {
      # Build interpolating model
      pred <- predict(object=model, newdata=model@X, type=type, checkNames = FALSE)
      mk <- pred$mean
      optim.param$ymin <- min(mk)
      if(exists("optim.model")) rm(optim.model)
      
      try(optim.model <- km(formula=model@trend.formula, design=model@X, response=mk,
    	            covtype=model@covariance@name, coef.trend=model@trend.coef, coef.cov=covparam2vect(model@covariance),
    	            coef.var=model@covariance@sd2,control=model@control))

      # If km crashes: add nugget to the interpolating model
      if(!exists("optim.model"))
      { print("Error occured during model update - small nugget added to the reinterpolating model")
        optim.model <- km(formula=model@trend.formula, design=model@X, response=mk,
    	            covtype=model@covariance@name, coef.trend=model@trend.coef, coef.cov=covparam2vect(model@covariance),
    	            coef.var=model@covariance@sd2, nugget=1e-8,control=model@control)
      }

      # Choose new observation based on EI
      oEGO <- max_EI(model=optim.model, plugin=min(optim.model@y), upper=upper, type=type, lower=lower, parinit=parinit, control=control)
      x.new <- oEGO$par
      i.best <- 0
    }
    #-------------------------------------------------------------------------------------------------------
    #-- Find the new observation based on the EI WITH PLUGIN criterion -------------------------------------
    #-------------------------------------------------------------------------------------------------------
    else if (optim.crit =="EI.plugin")
    {
      # Set plugin value (depending on "optim.param$plugin.type": "quantile", "ytilde" or "other")
      if (is.null(optim.param$plugin.type))
      { cat("Plugin type not provided: default value ytilde is used \n")
        optim.param$plugin.type <- "ytilde"}
      
      if (optim.param$plugin.type=="quantile")
      {  pred <- predict(object=model, newdata=model@X, type=type, checkNames = FALSE)
         mk <- pred$mean
         sk <- pred$sd
         if (is.null(optim.param$quantile))
         { cat("Quantile level not provided: default value 0.5 is used \n")
           optim.param$quantile <- 0.5
         }
         qk <- mk + qnorm(optim.param$quantile)*sk
         plugin <- min(qk)
      } else if (optim.param$plugin.type=="ytilde")
      {  plugin <- min(model@y)
      } else if (optim.param$plugin.type=="fixed")
      {  plugin <- optim.param$plugin
      } else 
      { cat("Unknown plugin type: default value ytilde is used \n")
        optim.param$plugin.type="ytilde"
        plugin <- min(model@y)
      }

      # Choose new observation based on EI.plugin
      oEGO <- max_EI(model=model, plugin=plugin, type=type, lower=lower, upper=upper, parinit=parinit, control=control)
      x.new <- oEGO$par
      EI <- oEGO$val

      # Compare with values at DoE
      EI.doe <- rep(0, 1, model@n)
      for (i in 1:model@n)
      {  EI.doe[i] <- EI(x=model@X[i,], model=model, type=type, plugin=plugin)}
      i.best <- which.max(EI.doe)

      # Change new observation if DoE is better
      if (EI.doe[i.best] > EI)
      {  x.new <- t(as.numeric(model@X[i.best,]))
         add.obs <- FALSE        }
    }
    #-------------------------------------------------------------------------------------------------------
    #-- Find the new observation based on the EXPECTED QUANTILE IMPROVEMENT criterion ----------------------
    #-------------------------------------------------------------------------------------------------------
    else if (optim.crit =="EQI")
    {
      # Set future noise variance, compute current minimal quantile
      if (is.null(optim.param))
      { beta <- 0.9
      } else { beta <- optim.param$quantile }

      new.noise.var <- noise.var/(n.ite+1-i.time.steps)
      pred <- predict(object=model, newdata=model@X, type=type, checkNames = FALSE)
      qk <- pred$mean + qnorm(beta)*pred$sd
      q.min <- min(qk)

      # Choose new observation based on EQI
      oEGO <- max_EQI(model=model, new.noise.var=new.noise.var, type=type,
                              beta=beta, q.min=q.min, lower=lower, upper=upper, parinit=parinit, control=control)
      x.new <- oEGO$par
      EQI.global.search <- oEGO$val

      # Compare with values at DoE
      EQI.doe <- rep(0, 1, model@n)
      for (i in 1:model@n)
      {  EQI.doe[i] <- EQI(x=model@X[i,], model=model, new.noise.var=new.noise.var, q.min=q.min, beta=beta, type=type)}
      i.best <- which.max(EQI.doe)

      # Change new observation if DoE is better
      if (EQI.doe[i.best] > EQI.global.search)
      {  x.new <- t(as.numeric(model@X[i.best,]))
         add.obs <- FALSE        }
    }
    #-------------------------------------------------------------------------------------------------------
    #-- Find the new observation based on the MINIMAL QUANTILE criterion -----------------------------------
    #-------------------------------------------------------------------------------------------------------
    else if (optim.crit =="min.quantile")
    {
      if (is.null(optim.param))
      { beta <- 0.1 
      } else { beta <- optim.param$quantile }

      # Choose new observation based on minimum quantile
      oEGO <- min_quantile(model=model, beta=beta, type=type, lower=lower, upper=upper, parinit=parinit, control=control)
      x.new <- oEGO$par
      q <- oEGO$val

      # Compare with values at DoE
      pred <- predict(object=model, newdata=model@X, type=type, checkNames = FALSE)
      mk <- pred$mean
      sk <- pred$sd
      q.doe <- mk + qnorm(beta)*sk
      i.best <- which.min(q.doe)

      # Change new observation if DoE is better
      if (q.doe[i.best] < q)
      {  x.new <- t(as.numeric(model@X[i.best,]))
         add.obs <- FALSE        }
    }
    #-------------------------------------------------------------------------------------------------------
    #-- Find the new observation based on the AUGMENTED EXPECTED IMPROVEMENT criterion ---------------------
    #-------------------------------------------------------------------------------------------------------
    else if (optim.crit =="AEI")
    {
      if (is.null(optim.param))
      { beta <- 0.75
      } else { beta <- optim.param$quantile }

      # Find current best and y.min
      pred <- predict(object=model, newdata=model@X, type=type, checkNames = FALSE)
      mk <- pred$mean
      sk <- pred$sd
      qk <- mk + qnorm(beta)*sk
      y.min <- mk[which.min(qk)]

      # Choose new observation based on AEI
      oEGO <- max_AEI(model=model, y.min=y.min, new.noise.var=noise.var, lower=lower, upper=upper, parinit=parinit, control=control, type=type)
      x.new <- oEGO$par
      AEI <- oEGO$val

      # Compare with values at DoE
      AEI.doe <- rep(0, 1, model@n)
      for (i in 1:model@n)
      {  AEI.doe[i] <- AEI(x=model@X[i,], model=model, y.min=y.min, type=type, new.noise.var=noise.var)}
      i.best <- which.max(AEI.doe)

      # Change new observation if DoE is better
      if (AEI.doe[i.best] > AEI)
      {  x.new <- t(as.numeric(model@X[i.best,]))
         add.obs <- FALSE        }
    }

    #-------------------------------------------------------------------------------------------------------
    #-- Find the new observation based on the APPROXIMATE KNOWLEDGE GRADIENT criterion ---------------------
    #-------------------------------------------------------------------------------------------------------
    else if (optim.crit =="AKG")
    {
      # Choose new observation based on AKG
      oEGO <- max_AKG(model=model, new.noise.var=noise.var, type=type, lower=lower, upper=upper, parinit=parinit, control=control)
      x.new <- oEGO$par
      AKG <- oEGO$val

      # Compare with values at DoE
      AKG.doe <- rep(0, 1, model@n)
      for (i in 1:model@n)
      {  AKG.doe[i] <- AKG(x=model@X[i,], model=model, type=type, new.noise.var=noise.var)}
      i.best <- which.max(AKG.doe)

      # Change new observation if DoE is better
      if (AKG.doe[i.best] > AKG)
      {  x.new <- t(as.numeric(model@X[i.best,]))
         add.obs <- FALSE        }
    }
    #-------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------
    #-- Generate observation, update model -----------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------
    # Generate observation
    y.new <- funnoise(x.new)

    # Update model
    if (add.obs)
    { cat(c("Creating obs:", as.numeric(x.new),"\n"))
    } else
    { cat(c("Repeating obs:", as.numeric(x.new),"\n"))}
 
    upmod <- update_km_noisyEGO(model=model, x.new=x.new, y.new=y.new, noise.var=noise.var, type=type,
                                add.obs=add.obs, index.in.DOE=i.best, CovReEstimate=CovReEstimate, NoiseReEstimate=NoiseReEstimate, 
                                estim.model=estim.model, nugget.LB=nugget.LB)
    model <- upmod$model
    if (NoiseReEstimate)
    { estim.model <- upmod$estim.model
      noise.var   <- upmod$noise.var
    }
    
    #-------------------------------------------------------------------------------------------------------
    #-- Save X path, observations and thetas ---------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------
    mu    <- model@trend.coef
    sd2   <- model@covariance@sd2
    range <- model@covariance@range.val
    theta <- c(mu, sd2, range)

    if (i.time.steps==1)
    {  all.X <- t(x.new)
       all.y <- t(y.new)
       all.thetas <- theta
       all.noise.var <- noise.var
    } else
    {  all.X <- cbind(all.X, t(x.new))
       all.y <- c(all.y, t(y.new))
       all.thetas <- cbind(all.thetas, theta)
       all.noise.var <- c(all.noise.var, t(noise.var))
    }
  }
  ### MAIN LOOP ENDS #######################################################################################
  ##########################################################################################################
  
  ############################################################################################################
  ######## COMPUTE BEST POINT ################################################################################
  ############################################################################################################
  # All other config
  
  if (!exists("optim.param$quantile")) { beta <- 0.5 
  } else { beta <- optim.param$quantile }
  
  pred <- predict(object=model, newdata=model@X, type=type, checkNames = FALSE)
  mk <- pred$mean
  sk <- pred$sd
  qk <- mk + qnorm(beta)*sk
  
  i.best <- which.min(qk)
  x.best <- model@X[i.best,]
  y.best <- model@y[i.best]
  
  if (optim.crit=="EI.plugin") 
  { if (optim.param$plugin.type=="ytilde")
    { # EI.plugin with ytilde treated differently
      i.best <- which.min(model@y)
      x.best <- model@X[i.best,]
      y.best <- model@y[i.best]
    }
  } 
  
  ############################################################################################################
  ######## RETURN OPTIMIZATION RESULTS #######################################################################
  ############################################################################################################
  optim.result$model <- model
  
  optim.result$best.x <- x.best
  optim.result$best.y <- y.best
  optim.result$best.index <- i.best
  
  optim.result$history.hyperparam <- all.thetas
  optim.result$history.x <- all.X
  optim.result$history.y <- all.y
  
  if (NoiseReEstimate)
  { optim.result$estim.model <- estim.model
    optim.result$history.noise.var <- all.noise.var}

  return(optim.result)
}