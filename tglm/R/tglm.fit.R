# library(mvtnorm);
# library(truncnorm);
# library(coda);
# library(BayesLogit);

is.binary = function(x){
  length(unique(x)) == 2;
}

##############################################################
######   Main function: logistic / probit regression    ######
##############################################################
######################## Inputs ########################
## y: vector of length n. Binary response.
## X: n * p matrix. Design matrix.
## iter: number of iterations in the MCMC.
## thin: thinning. Save one iteration for every "thin" number of iterations.
## burnin: the ratio of the burn-in length (over the total length).
## method: 'logistic' or 'probit'
## save.latent: logical. Whether to save the MCMC for the latent variable 
##             (omega for logistic regression, or z for probit regression). 
##             Since its length is n, it takes too much space when n is large.
## df: degree freedom of the independent t prior on both intercept and slopes. 
##     When df = Inf, use independent normal prior.
## slope.scale: the scale (or sd) parameter of the t (or normal) prior on slopes.
##				can be a scalar or a p-dim vector.
## intercept.scale: the scale (or sd) parameter of the t (or normal) prior 
##                  on the intercept.
## center.binary: logical. Whether to center binary predictors. If TRUE, the posterior
##                samples of beta are post-processed to the orginal scale.
## center.continuous: logical. Whether to center and rescale the non-binary predictors. 
##                    If TRUE, the posterior samples of beta are post-processed 
##                    to the orginal scale.

######################## Outputs ########################
## Inference: a matrix. Posterior mean, posterior sd, and 95% HPD interval of beta. 
##            (including both intercept and slopes).
## Beta: (iter / thin) * (p + 1) matrix. Posterior samples of beta.
## Gamma: (iter / thin) * (p + 1) matrix. Posterior samples of gamma (the latent 
##        parameters in the hiearchial preserntation of t priors). Included in the 
##        output only if df < Inf.
## Omega: (iter / thin) * n matrix. Posterior samples of omega (the latent 
##        parameters in the Polya Gamma hiearchial preserntation logistic likelihood). 
##        Included in the output only if save.omega = TRUE.
#########################################################

tglm.fit = function(y, X, iter = 1e5, thin = max(1, round(iter/2e3)), burnin = 0.5, 
                  method = 'logistic', df = 1, slope.scale = 2.5, intercept.scale = 10,
                  save.latent = FALSE, center.binary = TRUE, scale.continuous = TRUE, 
                  beta.original = TRUE, track.time = TRUE, show.summary = TRUE){

  t.start = proc.time();
  #### normalizing X ######
  n = dim(X)[1];
  p = dim(X)[2];
  Xmean = apply(X, 2, mean);
  Xsd = apply(X, 2, sd);
  Xbinary = apply(X, 2, is.binary);

  ## center binary predictors
  if(center.binary == TRUE)
  	X[, Xbinary] = sweep(as.matrix(X[, Xbinary]), 2, Xmean[Xbinary]);
  ## center and scale continuous predictors
  if(scale.continuous == TRUE){
  	X[, !Xbinary] = sweep(as.matrix(X[, !Xbinary]), 2, Xmean[!Xbinary]);
  	X[, !Xbinary] = sweep(as.matrix(X[, !Xbinary]), 2, Xsd[!Xbinary] * 2, FUN = "/");
  }
  
  ## include intercept column (all 1's) in the design matrix
  X = cbind(1, X);
  
  ## sigma: a (p + 1)-vector of prior scales
  if(length(slope.scale) == 1)
    sigma = c(intercept.scale, rep(slope.scale, p));
  if(length(slope.scale) == p)
    sigma = c(intercept.scale, slope.scale);
  if(length(slope.scale) != p && length(slope.scale) != 1)
    stop('Error: length of slope.scale should equal 1 or number of columns of X!');
    
  
  #### intial values of beta (including intercept), gamma (from their priors) ####
  beta = rnorm(p + 1, 0, sigma);
  if(df < Inf)
    gamma = 1 / rgamma(p + 1, df / 2, rate = sigma^2 * df / 2);
  if(df == Inf)
    gamma = sigma^2;
  if(method == 'logistic')
    omega = rpg(n, 1, X %*% beta);
  if(method == 'probit'){
    ## upper bounds and lower bounds of z (depend on y)
    z.lower = z.upper = rep(0, n);
    z.lower[y == 0] = -Inf;
    z.upper[y == 1] = Inf;
    z = rtruncnorm(n, a = z.lower, b = z.upper, mean = X %*% beta, sd = 1);  	
  }

  #### Save MCMC chains ####
  Beta = matrix(NA, nrow = iter / thin, ncol = p + 1);
  if(df < Inf)
    Gamma = matrix(NA, nrow = iter / thin, ncol = p + 1);
  if(save.latent == TRUE){
  	if(method == 'logistic')
      Omega = matrix(NA, nrow = iter / thin, ncol = n);
    if(method == 'probit')
      Z = matrix(NA, nrow = iter / thin, ncol = n);
  }

  #### Gibbs sampler ####
  for(it in 1:iter){

    if(method == 'logistic'){    
      ## update beta ##
      cov.beta = solve(t(X) %*% sweep(X, 1, omega, FUN = "*") + diag(1 / gamma));
      beta = c(rmvnorm(1, cov.beta %*% (t(X) %*% (y - 0.5)), cov.beta));
      ## update omega ##
      omega = rpg(n, 1, X %*% beta);
    }
    
    if(method == 'probit'){    
      ## update beta ##
      cov.beta = solve(t(X) %*% X + diag(1 / gamma));
      beta = c(rmvnorm(1, cov.beta %*% (t(X) %*% z), cov.beta));
      ## update z ##
      z = rtruncnorm(n, a = z.lower, b = z.upper, mean = X %*% beta, sd = 1);
    }
     
    ## update gamma ##
    if(df < Inf)
      gamma = 1 / rgamma(p + 1, (df + 1) / 2, rate = (beta^2 + sigma^2 * df) / 2);
        
    ## save every thin iterations
    if(it %% thin == 0){
  	  record = it / thin;
  	  Beta[record, ] = beta;
  	  if(df < Inf)
  	    Gamma[record, ] = gamma;  	    
      if(save.latent == TRUE){
  	    if(method == 'logistic')
          Omega[record, ] = omega;
        if(method == 'probit')
          Z[record, ] = z;
      }
      
      
      ## show: x0% completed
      if( (it * 10) %% iter == 0 && it != iter && track.time == TRUE)
        cat(paste( (it * 100) / iter), '% completed...\n', sep = '');
      if( it == iter && track.time == TRUE){
        cat(paste( (it * 100) / iter), '% completed.\n', sep = '');
        cat('\n');
      } 

    }
    
  }
  
  ## track time
  if(track.time == TRUE){
    t.finish = proc.time();
    cat('Time used (in second): \n')
    print(t.finish - t.start);
    cat('\n');
  }
  
  ## transform beta back to orginal scale (before centering and scaling)
  if(center.binary == TRUE && beta.original == TRUE)
    Beta[, 1] = Beta[, 1] - as.matrix(Beta[, which(Xbinary) + 1]) %*% Xmean[Xbinary];
  if(scale.continuous == TRUE && beta.original == TRUE){
  	Beta[, 1] = Beta[, 1] - as.matrix(Beta[, which(!Xbinary) + 1]) %*% (Xmean[!Xbinary] / 2 / Xsd[!Xbinary]);
  	Beta[, which(!Xbinary) + 1] = as.matrix(Beta[, which(!Xbinary) + 1]) / 2 / Xsd[!Xbinary];
  }

  ##########################################
  ######  Posterior Inference:        ######
  ######  mean, sd, HDP               ######
  ##########################################
  
  keep = round(iter / thin * burnin): (iter / thin);
  varnames = colnames(X[, -1]);
  if(length(varnames) == 0)
    varnames = paste('X', 1:p, sep = '');
    
  inference = matrix(NA, ncol = 4, nrow = p + 1);
  rownames(inference) = c('intercept', varnames);
  colnames(inference) = c('Est', 'Std', '95HPDlower','95HPDupper');
  
  ## posterior mean and HPD
  inference[, 1] = apply(as.matrix(Beta[keep, ]), 2, mean);
  inference[, 2] = apply(as.matrix(Beta[keep, ]), 2, sd);
  inference[, 3:4] = HPDinterval(as.mcmc(Beta[keep, ]), 0.95);
    
  ## show summary: inference
  if(show.summary == TRUE){
    cat('Posterior inference: \n');
    print(round(inference, 4));
  }

   
  if(save.latent == FALSE && df < Inf)
    return(list(Beta = Beta, Gamma = Gamma, inference = inference));
  if(save.latent == FALSE && df == Inf)
    return(list(Beta = Beta, inference = inference));
  if(save.latent == TRUE && df < Inf){
  	if(method == 'logistic')
      return(list(Beta = Beta, Gamma = Gamma, Z = Omega, inference = inference));
  	if(method == 'probit')
      return(list(Beta = Beta, Gamma = Gamma, Z = Z, inference = inference));
  }
  if(save.latent == TRUE && df == Inf){
  	if(method == 'logistic')
      return(list(Beta = Beta, Z = Omega, inference = inference));
  	if(method == 'probit')
      return(list(Beta = Beta, Z = Z, inference = inference));
  }
 
}


