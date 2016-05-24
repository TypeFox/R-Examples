gaussian.CARadaptive <- function(formula, data = NULL, W, burnin, n.sample, thin = 1, prior.mean.beta = NULL, prior.var.beta = NULL, prior.tau2 = NULL, prior.nu2 = NULL, rhofix = NULL, epsilon = 0, verbose = TRUE)
    { 
  #### Check on the verbose option
  if(is.null(verbose))     verbose = TRUE     
  if(!is.logical(verbose)) stop("the verbose option is not logical.", call.=FALSE)
  if(verbose) cat("Setting up the model\n"); a<-proc.time()
  
  blocksize.beta <- 5
  blocksize.v    <- 10
  z              <- which(W > 0, arr.ind = T)
  locs           <- z[which(z[,1] < z[,2]), ]
  char.locs      <- paste(locs[,1], ".", locs[,2], sep = "")
  n.edges        <- nrow(locs)
  
  # convert the supplied adjacency matrix into a spam matrix, if required.
  if(class(W) == "matrix") W <- as.spam(W)
  if(!class(W) %in% c("matrix", "spam")) stop("W must be an object with class \"matrix\" or \"spam\"", call.=FALSE)  
  
  logit     <- function(p) log(p/(1-p))
  inv_logit <- function(v) 1/(1+exp(-v))
  
  # interpret the formula
  frame <- try(suppressWarnings(model.frame(formula, data = data, na.action=na.pass)), silent=TRUE)
  if(class(frame)=="try-error") stop("the formula inputted contains an error, e.g the variables may be different lengths.", call.=FALSE)
  X <- try(suppressWarnings(model.matrix(object=attr(frame, "terms"), data=frame)), silent=TRUE)
  if(class(X)=="try-error") stop("the covariate matrix contains inappropriate values.", call.=FALSE)
  if(sum(is.na(X))>0) stop("the covariate matrix contains missing 'NA' values.", call.=FALSE)
  # get summaries of the model matrix
  p       <- ncol(X)
  y       <- model.response(frame)
  n.sites <- as.integer(nrow(W))
  n.time  <- as.integer(length(y)/n.sites)
  k       <- as.integer(round(n.sites*n.time, 0))
  
  # check y, the model response
  if(sum(is.na(y))>0) stop("the response has missing 'NA' values.", call.=FALSE)
  if(!is.numeric(y)) stop("the response variable has non-numeric values.", call.=FALSE)
  
  #### Check and specify the priors
  if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
  if(length(prior.mean.beta)!=p) stop("the vector of prior means for beta is the wrong length.", call.=FALSE)    
  if(!is.numeric(prior.mean.beta)) stop("the vector of prior means for beta is not numeric.", call.=FALSE)    
  if(sum(is.na(prior.mean.beta))!=0) stop("the vector of prior means for beta has missing values.", call.=FALSE)       
  
  if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
  if(length(prior.var.beta)!=p) stop("the vector of prior variances for beta is the wrong length.", call.=FALSE)    
  if(!is.numeric(prior.var.beta)) stop("the vector of prior variances for beta is not numeric.", call.=FALSE)    
  if(sum(is.na(prior.var.beta))!=0) stop("the vector of prior variances for beta has missing values.", call.=FALSE)    
  if(min(prior.var.beta) <=0) stop("the vector of prior variances has elements less than zero", call.=FALSE)
  prior.prec.beta <- diag.spam(1/prior.var.beta, nrow = length(prior.var.beta))
  
  if(is.null(prior.tau2)) prior.tau2 <- c(0.001, 0.001)
  if(length(prior.tau2)!=2) stop("the prior value for tau2 is the wrong length.", call.=FALSE)    
  if(!is.numeric(prior.tau2)) stop("the prior value for tau2 is not numeric.", call.=FALSE)    
  if(sum(is.na(prior.tau2))!=0) stop("the prior value for tau2 has missing values.", call.=FALSE)    
  
  if(is.null(prior.nu2)) prior.nu2 <- c(0.001, 0.001)
  if(length(prior.nu2)!=2) stop("the prior value for nu2 is the wrong length.", call.=FALSE)    
  if(!is.numeric(prior.nu2)) stop("the prior value for nu2 is not numeric.", call.=FALSE)    
  if(sum(is.na(prior.nu2))!=0) stop("the prior value for nu2 has missing values.", call.=FALSE)   
  
  # identify and error check the offset term, if it exists.
  offset <- try(model.offset(frame), silent=TRUE)
  if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
  if(is.null(offset))              offset <- rep(0,(n.time * n.sites))
  if(sum(is.na(offset))>0)         stop("the offset has missing 'NA' values.", call.=FALSE)
  if(!is.numeric(offset))          stop("the offset variable has non-numeric values.", call.=FALSE) 
  
  #### Format and check the MCMC quantities
  if(is.null(burnin)) stop("the burnin argument is missing", call.=FALSE)
  if(is.null(n.sample)) stop("the n.sample argument is missing", call.=FALSE)
  if(!is.numeric(burnin)) stop("burn-in is not a number", call.=FALSE)
  if(!is.numeric(n.sample)) stop("n.sample is not a number", call.=FALSE) 
  if(!is.numeric(thin)) stop("thin is not a number", call.=FALSE)
  if(n.sample <= 0) stop("n.sample is less than or equal to zero.", call.=FALSE)
  if(burnin < 0) stop("burn-in is less than zero.", call.=FALSE)
  if(thin <= 0) stop("thin is less than or equal to zero.", call.=FALSE)
  if(n.sample <= burnin)  stop("Burn-in is greater than n.sample.", call.=FALSE)
  if(n.sample <= thin)  stop("thin is greater than n.sample.", call.=FALSE)
  if(burnin!=round(burnin)) stop("burnin is not an integer.", call.=FALSE) 
  if(n.sample!=round(n.sample)) stop("n.sample is not an integer.", call.=FALSE) 
  if(thin!=round(thin)) stop("thin is not an integer.", call.=FALSE) 
  
  ## Check for linearly related columns
  cor.X <- suppressWarnings(cor(X))
  diag(cor.X) <- 0
  
  if(max(cor.X, na.rm=TRUE)==1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
  if(min(cor.X, na.rm=TRUE)==-1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
  
  if(p>1)
  {
    if(sort(apply(X, 2, sd))[2]==0) stop("the covariate matrix has two intercept terms.", call.=FALSE)
  }else
  {
  }
  

  
  ## Standardise the model matrix, 
  X.standardised   <- X
  X.sd             <- apply(X, 2, sd)
  X.mean           <- apply(X, 2, mean)
  X.indicator      <- rep(NA, p)       # To determine which parameter estimates to transform back
  for(j in 1:p){
    if(length(table(X[ ,j])) > 2){
      X.indicator[j] <- 1
      X.standardised[ ,j] <- (X[ ,j] - mean(X[ ,j])) / sd(X[ ,j])
    }else if(length(table(X[ ,j]))==1){
      X.indicator[j] <- 2
    }else{
      X.indicator[j] <- 0
    }
  } 
  
  # based on the blocksize.v provided create lists with relevent bits for untransformed edge parameter update
  if(is.numeric(blocksize.v)){
    ## Compute the blocking structure for v
    fromto     <- seq(0, n.edges, by = blocksize.v)
    fromto[1]  <- 0
    if(!n.edges %in% fromto) fromto <- c(fromto, n.edges)
    n.blocks   <- length(fromto) - 1
    blockinds  <- vector("list", length = n.blocks)
    for(i in 1:n.blocks)  blockinds[[i]]    <- (fromto[i] + 1):fromto[i + 1]
  } 
  
  # propose starting values for the adjacency elements (very close to 1)
  # current temporary version of the adacency is W_current
  v                                <- logit(rtrunc(n.edges, spec = "norm", mean = 0.999, sd = 0.001, a = 0, b=1))
  v_15                             <- v - 15
  vqform_current                   <- sum(v_15^2)
  W_current                        <- W
  W_current[locs][1:n.edges]       <- inv_logit(v)
  W_current[locs[,2:1]][1:n.edges] <- inv_logit(v)
  
  # given the temporary adjacency, construct a temporary (Q.space) and proposal (Q.space.prop)
  # for the prior ICAR precision for phi.  Associated with these is the triplet form tripList.
  # get the cholesky of Q.space, and its determinant
  # if rho is not fixed, then ridge must be fixed
  rho                     <- ifelse(!is.null(rhofix), rhofix, 0.99)
  fixedridge              <- epsilon
  if(rho==1) fixedridge <- 0.0001
  tripList                <- vector("list", length = 2)
  tripList[[1]]           <- cbind(1:nrow(W_current), 1:nrow(W_current), rowSums(W_current) + fixedridge)
  tripList[[2]]           <- cbind(rbind(locs, locs[,2:1]), -rep(inv_logit(v), 2))
  Q.space.trip            <- rbind(tripList[[1]], tripList[[2]])
  Q.space.trip            <- updatetriplets_rho(trips = Q.space.trip, nsites = n.sites, rho_old = 1, rho_new = rho, fixedridge = fixedridge)  
  Q.space <- Q.space.prop <- spam(list(i = Q.space.trip[,1], j = Q.space.trip[,2], Q.space.trip[,3]))
  chol.Q.space            <- chol.spam(Q.space)
  Q.space.det.old         <- n.time*2*determinant(chol.Q.space, logarithm = T)$modulus
  
  # propose an initial value for alpha, the temporal correlation parameter
  # using alpha, create initial temporal precision matrices Q.time
  alpha <- 1
  if(n.time > 1){
    # this bit constructs Q.time, temporal precision, its determinant, and triplet form
    Q.block      <- as.spam(crossprod(diff(diag(n.time))))
    Q.block[1,1] <- Q.block[1,1] + 1
    Dg           <- diag.spam(diag.spam(Q.block))
    R            <- Q.block - Dg
    Dtime        <- diag.spam( c(rep(1,nrow(Q.block)-1), 0))
    Dg           <- Dg - Dtime
    Q.time       <- Dg + Dtime*alpha^2+ R*alpha
    Q.time[n.time,n.time] <- 1
    Q.det        <- determinant(Q.time, logarithm = T)
    detTime      <- as.numeric(0.5*n.sites*(Q.det$m)*(Q.det$s))
    Q.time.trip  <- Reduce("cbind", triplet(Q.time)) 
  }  else {
    # if n.time == 1, detTime equals 1 and Q.time is just a 1 times 1 matrix.
    Q.time       <- 1
    detTime      <- 1
    Q.time.trip  <- matrix(rep(1, 3), ncol = 3) 
  }
  
  # MCMC parameter starting values
  W.tune        <- 1
  rho.tune      <- 0.1
  tau_v         <- 200
  prior.max.tau <- 1000
  increment     <- 0
   mod.glm <- glm(y~X.standardised-1, offset=offset)
  beta.mean <- mod.glm$coefficients
  beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
  beta_par <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
  res.temp <- y - X.standardised %*% beta_par - offset
  res.sd <- sd(res.temp, na.rm=TRUE)/5
  phi <- rnorm(n=k, mean=0, sd = res.sd)
  tau <- var(phi)/10
  nu2 <- var(phi)/10
  phiQphi       <- qform_ST(Qspace = Q.space.trip, Qtime = Q.time.trip, phi = phi, nsites = n.sites) 
  XB            <- X.standardised %*% beta_par
  tau_v.shape   <- (n.edges/2) +  prior.tau2[1]
  tau_phi_shape <- (n.sites*n.time/2) + prior.tau2[1]
  # general MCMC housekeeping
  n.save        <- ifelse(thin == 1, (n.sample - burnin), (n.sample - burnin) / thin)
  accept.all    <- rep(0, 8)
  accept        <- accept.all
  # storage of parameters in the MCMC, 
  samples.beta  <- array(NA, c(n.save, p))
  samples.phi   <- array(NA, c(n.save, n.sites * n.time))
  samples.tau2  <- samples.deviance <- samples.vtau2 <- samples.alpha <- samples.nu2 <- samples.rho <- matrix(0, n.save, 1)
  samples.v     <- matrix(0, ncol = n.edges, nrow = c(n.save, n.sites*n.time))
  samples.fit   <- array(NA, c(n.save, n.sites * n.time))  
  samples.like <- array(NA, c(n.save, n.sites*n.time))  
  
  # turn off spam check options to speed things up (a bit)
  spam.options( "cholsymmetrycheck" = FALSE)
  spam.options( "cholpivotcheck" = FALSE)
  spam.options( "safemode" = c(F, F, F))
  
  ## Compute the blocking structure for beta     
  if(blocksize.beta >= p){
    n.beta.block <- 1
    beta.beg <- 1
    beta.fin <- p
  } else {
    n.standard <- 1 + floor((p-blocksize.beta) / blocksize.beta)
    remainder <- p - n.standard * blocksize.beta 
    if(remainder==0){
      beta.beg <- c(1,seq((blocksize.beta+1), p, blocksize.beta))
      beta.fin <- seq(blocksize.beta, p, blocksize.beta)
      n.beta.block <- length(beta.beg)
    } else {
      beta.beg <- c(1, seq((blocksize.beta+1), p, blocksize.beta))
      beta.fin <- c(seq((blocksize.beta), p, blocksize.beta), p)
      n.beta.block <- length(beta.beg)
    }
  }         
  corr.beta       <- t(X.standardised) %*% X.standardised
  corr.beta.inv   <- solve(t(X.standardised) %*% X.standardised)
  chol.corr.beta  <- chol(corr.beta.inv)    
  
  # the perm ordering is used to map the @entries slot ordering to the ordering used when 'triplet' is called
  perm           <- order(Q.space.trip[,1], Q.space.trip[,2])
  diags.space    <- which(Q.space.trip[perm,1] == Q.space.trip[perm,2])
  if(n.time > 1) diag.time      <- Reduce("cbind", triplet(diag.spam(n.time - 1)))
  time.last.diag <- which((Q.time.trip[,1] == Q.time.trip[,2]) & (Q.time.trip[,1] == n.time))
  lastblock      <- (k - n.sites+1):k
  firstblock     <- 1:n.sites
  
  ## Start timer
  if(verbose){
    cat("Generating", n.sample, "samples\n", sep = " ")
    progressBar            <- txtProgressBar(style = 3)
    percentage.points      <- round((1:100/100)*n.sample)
  } else percentage.points <- round((1:100/100)*n.sample)     
  
  
  #   -------------------------------------------------------------------------------------------
  #   START THE MCMC SAMPLING
  #   -------------------------------------------------------------------------------------------
  
  for(j in 1:n.sample){ 
    # START ITERATING, ONLY SAVE thin^th ITERATION
    save.iter <- j > burnin && ((j %% thin == 0) | thin == 0)
    if(save.iter) increment <- increment+1
    
    
    # update ALPHA
    if(n.time > 1){
      phifirst         <- phi[-firstblock]
      philast          <- phi[-lastblock]
      philastQphilast  <- qform_ST(Qspace = Q.space.trip, Qtime = diag.time, phi = philast, nsites = n.sites)   
      phifirstQphilast <- qform_ST_asym(Qspace = Q.space.trip, Qtime = diag.time, phi1 = phifirst, phi2 = philast, nsites = n.sites) 
      mu_alpha         <- phifirstQphilast/philastQphilast
      mu_sigmasq       <- tau/philastQphilast
      alpha            <- rtrunc(n=1, spec="norm", a=10^-5, b=1 - 10^-5,  mean=mu_alpha, sd = sqrt(mu_sigmasq))
      Q.time.trip      <- update_Qtime(Q.time.trip, alpha, time.last.diag - 1)
      phiQphi          <- qform_ST(Qspace = Q.space.trip, Qtime = Q.time.trip, phi = phi, nsites = n.sites)   
      detTime          <- determinant(Q.time, logarithm = TRUE)
      detTime          <- (detTime$m)*(detTime$s)
    }
    
    # Gibbs update of tau_v
    tau_scale  <- vqform_current/2 + prior.tau2[2]
    tau_v      <- 1/rtrunc(n=1, spec="gamma", a=0.000001, b=Inf, shape=tau_v.shape, scale=(1/tau_scale))
    v.proposal <- rtrunc(n = n.edges, spec="norm", a=-15, b=15,  mean = v, sd = W.tune)
    for(i in 1:n.blocks){
      # propose new v for the i^th block
      vnew                                 <- v
      block_inds                           <- blockinds[[i]]
      vnew[block_inds]                     <- v.proposal[block_inds] 
      # update the spatial precision matrix using c++ loop.  
      # This is efficient because changes are only made where vnew differs from v
      # combine the result back into triplet matrix (Q.space.trip.prop), and spam matrix (Q.space.prop)
      tripUpdate                           <- updatetripList2(Q.space.trip, vold = v, vnew = vnew, nedges = n.edges, 
                                                              nsites = n.sites, block = block_inds, 
                                                              block_length = length(block_inds), fixedridge = fixedridge, rho = rho) 
      Q.space.trip.prop                    <- tripUpdate[[1]]
      Q.space.trip.diff                    <- tripUpdate[[2]]
      # combine the result back into triplet matrix (Q.space.trip.prop), and spam matrix (Q.space.prop)      
      Q.space.prop@entries                 <- Q.space.trip.prop[perm,3]
      # acceptance ratio requires calculation of phi'Q_prop phi - phi'Q phi.
      # do this quickly by taking the difference between old and new triplets and working out the 
      # difference directly.  Much faster than working out quadratic forms seperately.
      Q.space.trip.diff[, 3]<- Q.space.trip[, 3] - Q.space.trip.prop[,3]
      phiQphi_phiQphiNew   <- qform_difference_ST(Qtrip = Q.space.trip.diff, Qtime = Q.time.trip, phi = phi, nsites = n.sites)
      # update the cholesky of the precision matrix & calculate the determinant
      chol.Q.space.prop    <- update(chol.Q.space, x = Q.space.prop) 
      detSpace             <- 2*determinant(chol.Q.space.prop, logarithm = T)$modulus
      Q.space.det.prop     <- n.sites*detTime + n.time*detSpace
      v_15_prop            <- vnew - 15
      vqform_prop          <- sum(v_15_prop^2)
      acceptance           <- exp(0.5*(Q.space.det.prop - Q.space.det.old) + (1/(2*tau))*(phiQphi_phiQphiNew) 
                                  + 0.5*(1/tau_v)*(vqform_current - vqform_prop))
      accept[8]      <- accept[8] + (1/n.blocks)
      if(runif(1)  <= acceptance){
        vqform_current   <- vqform_prop
        v                <- vnew
        accept[7]        <- accept[7] + (1/n.blocks)
        Q.space.det.old  <- Q.space.det.prop
        Q.space.trip     <- Q.space.trip.prop
        chol.Q.space     <- chol.Q.space.prop
        Q.space          <- Q.space.prop
      }
    }
    
    # Gibbs update beta
    var_mat_full_con   <- solve(((corr.beta/nu2) + prior.prec.beta))
    mean_full_con      <- ((((y - phi - offset)/nu2) %*% X.standardised) + (prior.mean.beta/prior.var.beta)) %*% var_mat_full_con
    beta_par           <- mvrnorm(n = 1, mu = mean_full_con, Sigma = var_mat_full_con)
    XB                 <- X.standardised %*% beta_par
    
    # Gibbs update of the variance of the normal likelihood
    lik_comp <- sum( (y - XB - phi - offset)^2 )/2
    nu2      <- 1 / rgamma(1, (prior.nu2[1] + 0.5*k), scale=(1/(prior.nu2[2] + lik_comp)))     
    
    # update PHI using one at a time M-H sampling
    nneighbours   <- diag.spam(Q.space)
    W_current     <- diag(nneighbours) - as.matrix(Q.space)
    phi_update    <- SPTICARphiGaussian(W = W_current, nsites = n.sites, ntimes = n.time, 
                                        phi = phi, nneighbours = nneighbours, 
                                        tau = tau, y = as.vector(y - offset), 
                                        alpha = alpha, XB = XB, lik_var = nu2)    
    phi       <- phi_update[[1]]
    phi       <- phi - mean(phi)
    
    # update rho, the spatial leroux parameter
    if(!is.null(rhofix)){
      proposal.rho <- rhofix
    } else {
      proposal.rho           <- rtrunc(n = 1, spec="norm", a=0, b=1, mean = rho, sd = rho.tune) 
    }
    Q.space.trip.prop      <- updatetriplets_rho(trips = Q.space.trip, nsites = n.sites, rho_old = rho, rho_new = proposal.rho, fixedridge = fixedridge)   
    Q.space.prop@entries   <- Q.space.trip.prop[perm,3]
    Q.space.trip.diff[, 3] <- Q.space.trip[, 3] - Q.space.trip.prop[,3]
    phiQphi_phiQphiNew     <- qform_difference_ST(Qtrip = Q.space.trip.diff, Qtime = Q.time.trip, phi = phi, nsites = n.sites)
    # update the cholesky of the precision matrix & calculate the determinant
    chol.Q.space.prop      <- update(chol.Q.space, x = Q.space.prop) 
    detSpace               <- 2*determinant(chol.Q.space.prop, logarithm = T)$modulus
    Q.space.det.prop       <- n.sites*detTime + n.time*detSpace
    acceptance             <- exp(0.5*(Q.space.det.prop - Q.space.det.old) + (1/(2*tau))*(phiQphi_phiQphiNew))
    accept[6]              <- accept[6] + 1
    if(runif(1)  <= acceptance){
      accept[5]        <- accept[5] + 1
      Q.space.det.old  <- Q.space.det.prop
      Q.space.trip     <- Q.space.trip.prop
      chol.Q.space     <- chol.Q.space.prop
      Q.space          <- Q.space.prop
      rho              <- proposal.rho
    }
    
    
    # Gibbs update TAU using the gamma distribution
    phiQphi    <- qform_ST(Qspace = Q.space.trip, Qtime = Q.time.trip, phi = phi, nsites = n.sites)     
    tau_scale  <- phiQphi/2 + prior.tau2[2]
    tau        <- 1/rtrunc(n=1, spec="gamma", a=0.000001, b=Inf, shape=tau_phi_shape, scale=(1/tau_scale))
    
    # calculate the deviance
    fitted   <- as.numeric(X.standardised %*% beta_par) + phi + offset
    deviance.all      <- dnorm(y, mean = fitted, sd = rep(sqrt(nu2), k), log=TRUE)
    like <- exp(deviance.all)
    dev          <- -2 * sum(deviance.all)   
    
    
    
    # save samples if past burnin 
    if(save.iter){
      samples.beta[increment,]      <- beta_par
      samples.phi[increment,]       <- phi
      samples.fit[increment, ]      <- fitted
      samples.tau2[increment,]      <- tau
      samples.deviance[increment,]  <- dev
      samples.vtau2[increment,]     <- tau_v
      samples.v[increment,]         <- v
      samples.alpha[increment,]     <- alpha
      samples.rho[increment,]       <- rho
      samples.nu2[increment,]       <- nu2
      samples.like[increment, ] <- like
    }
    
    
        # adjust the acceptance rate if required
        if(j %% 100 == 0){
             accept.w <- 100 * accept[7] / accept[8]
            if(is.null(rhofix))
            {
                accept.rho <- 100 * accept[5] / accept[6]     
            }else
            {
                accept.rho <- 45 
            }
            
            #### w tuning parameter
            if(accept.w > 40)
            {
                W.tune <- W.tune + 0.1 * W.tune
            }else if(accept.w < 20)              
            {
                W.tune <- W.tune - 0.1 * W.tune
            }else
            {
            }   
            
            #### rho tuning parameter
            if(accept.rho > 50)
            {
                rho.tune <- min(rho.tune + 0.1 * rho.tune, 0.5)
            }else if(accept.rho < 40)              
            {
                rho.tune <- rho.tune - 0.1 * rho.tune
            }else
            {
            }  
            accept.all         <- accept.all + accept
            accept             <- accept*0
            
    }
    
    
    
    # print progress to the console
    if(j %in% percentage.points & verbose) setTxtProgressBar(progressBar, j/n.sample)
  }
  
  # end timer
  if(verbose) cat("\nSummarising results"); close(progressBar)  
  
  ###################################
  #### Summarise and save the results 
  ###################################
  ## Compute the acceptance rates
  accept.beta  <- 100
  accept.phi   <- 100
  accept.rho <- 100 * accept.all[5] / accept.all[6]
  accept.w     <- 100 * accept.all[7] / accept.all[8]
  accept.alpha <- 100
  if(!is.null(rhofix))
  {
      accept.final <- c(accept.beta, accept.phi, accept.w)
      names(accept.final) <- c("beta", "phi", "w")  
  }else
  {
      accept.final <- c(accept.beta, accept.phi, accept.rho,accept.w)
      names(accept.final) <- c("beta", "phi", "rho", "w")  
  }  
  
  # ## Compute information criterion (DIC, DIC3, WAIC)
  median.beta        <- apply(samples.beta, 2, median)
  regression.mat     <- matrix(X.standardised %*% median.beta, nrow = n.sites, ncol = n.time, byrow=FALSE)   
  median.phi         <- matrix(apply(samples.phi, 2, median), nrow = n.sites, ncol = n.time)
  median.nu2         <- median(samples.nu2)
  offset.mat         <- matrix(offset, nrow = n.sites, ncol = n.time, byrow=FALSE) 
  fitted.median      <- as.numeric(median.phi + regression.mat + offset.mat)
  deviance.fitted    <- -2 * sum(dnorm(y, mean = fitted.median, sd = rep(sqrt(median.nu2), k), log = TRUE))
  p.d <- median(samples.deviance) - deviance.fitted
  DIC <- 2 * median(samples.deviance) - deviance.fitted     
  
  
  #### Watanabe-Akaike Information Criterion (WAIC)
  LPPD <- sum(log(apply(samples.like,2,mean)), na.rm=TRUE)
  p.w <- sum(apply(log(samples.like),2,var), na.rm=TRUE)
  WAIC <- -2 * (LPPD - p.w)
  
  
  ## Compute the LMPL
  CPO <- rep(NA, (n.sites * n.time))
  for(j in 1:(n.sites * n.time)){
    CPO[j] <- 1/median((1 / dnorm(y[j], mean=samples.fit[ ,j], sd=sqrt(samples.nu2))))   
  }
  LMPL <- sum(log(CPO))    
  
  
  ## Create the Fitted values
  fitted.values      <- apply(samples.fit, 2, median)
  residuals          <- as.numeric(y) - fitted.values
  
  #### transform the parameters back to the origianl covariate scale.
  samples.beta.orig <- samples.beta
  for(r in 1:p)
  {
    if(X.indicator[r]==1)
    {
      samples.beta.orig[ ,r] <- samples.beta[ ,r] / X.sd[r]
    }else if(X.indicator[r]==2 & p>1)
    {
      X.transformed <- which(X.indicator==1)
      samples.temp <- as.matrix(samples.beta[ ,X.transformed])
      for(s in 1:length(X.transformed))
      {
        samples.temp[ ,s] <- samples.temp[ ,s] * X.mean[X.transformed[s]]  / X.sd[X.transformed[s]]
      }
      intercept.adjustment <- apply(samples.temp, 1,sum) 
      samples.beta.orig[ ,r] <- samples.beta[ ,r] - intercept.adjustment
    }
  }
  
  #### Create a summary object  
  samples.beta.orig       <- mcmc(samples.beta.orig)
  summary.beta            <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
  summary.beta            <- cbind(summary.beta, rep(n.save, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
  rownames(summary.beta)  <- colnames(X)
  colnames(summary.beta)  <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")

  
  summary.hyper           <- array(NA, c(5, 7))   
  summary.hyper[1,1:3]    <- quantile(samples.nu2, c(0.5, 0.025, 0.975))
  summary.hyper[2,1:3]    <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
  summary.hyper[3,1:3]    <- quantile(samples.rho, c(0.5, 0.025, 0.975))
  summary.hyper[4,1:3]    <- quantile(samples.alpha, c(0.5, 0.025, 0.975))
  summary.hyper[5,1:3]    <- quantile(samples.vtau2, c(0.5, 0.025, 0.975))
  rownames(summary.hyper) <- c("nu2", "tau2", "rho.S",  "rho.T", "tau2.w")     
  summary.hyper[1, 4:7]   <- c(n.save, 100, effectiveSize(mcmc(samples.nu2)), geweke.diag(mcmc(samples.nu2))$z)     
  summary.hyper[2, 4:7]   <- c(n.save, 100, effectiveSize(mcmc(samples.tau2)), geweke.diag(mcmc(samples.tau2))$z)    
  summary.hyper[3, 4:7]   <- c(n.save, accept.rho, effectiveSize(mcmc(samples.rho)), geweke.diag(mcmc(samples.rho))$z)   
  summary.hyper[4, 4:7]   <- c(n.save, accept.alpha, effectiveSize(mcmc(samples.alpha)), geweke.diag(mcmc(samples.alpha))$z)   
  summary.hyper[5, 4:7]   <- c(n.save, 100, effectiveSize(mcmc(samples.vtau2)), geweke.diag(mcmc(samples.vtau2))$z)   
  
  if(!is.null(rhofix))
  {
      summary.hyper[3, ] <- c(rep(rhofix, 3),rep(NA, 4))    
  }
  
  summary.results         <- rbind(summary.beta, summary.hyper)
  summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
  summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
  
  # convert v back to w, summarise and create a 'fitted' adjacency matrix
  samples.w <- inv_logit(samples.v)
  colnames(samples.w) <- char.locs
  get_prop_thresh <- function(v, thresh) as.numeric(!((sum(v < thresh)/length(v)) < 0.99))
  bdry99          <- apply(samples.w, 2, get_prop_thresh, thresh = 0.5)
  bdryMN          <- apply(samples.w, 2, mean)
  Wmn <- W99      <- matrix(NA, nrow = n.sites, ncol = n.sites)
  W99[locs]       <- bdry99
  Wmn[locs]       <- bdryMN
  W99[locs]       <- bdry99
  W99[locs[ ,c(2,1)]] <- bdry99
  Wmn[locs]       <- bdryMN
  Wmn[locs[ ,c(2,1)]] <- bdryMN  
  
  
  ## Compile and return the results
  modelfit <- c(DIC, p.d, WAIC, p.w, LMPL)
  names(modelfit) <- c("DIC", "p.d", "WAIC", "p.w", "LMPL")
  model.string    <- c("Likelihood model - Gaussian (identity link function)", 
                       "\nLatent structure model - Adaptive autoregressive CAR model\n")
  samples.tau2all <- cbind(samples.tau2, samples.vtau2)
  colnames(samples.tau2all) <- c("tau2", "tau2.w")
  if(is.null(rhofix))
  {
      samples.rhoext <- cbind(samples.rho, samples.alpha)
      colnames(samples.rhoext) <- c("rho.S", "rho.T")
  }else
  {
      samples.rhoext <- cbind(samples.alpha)
      names(samples.rhoext) <- c("rho.T")
  }
  
  samples         <- list(beta = mcmc(samples.beta.orig), phi = mcmc(samples.phi), rho = mcmc(samples.rhoext), 
                          tau2 = mcmc(samples.tau2all), nu2 = mcmc(samples.nu2),  w = mcmc(samples.w), fitted = mcmc(samples.fit))  
  localised.structure <- list(Wmedian = Wmn, W99 = W99)
  results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=localised.structure, formula=formula, model=model.string,  X=X)
  class(results) <- "carbayesST"
  if(verbose)
  {
    b<-proc.time()
    cat(" finished in ", round(b[3]-a[3], 1), "seconds")
  }else
  {}
  return(results)
}







