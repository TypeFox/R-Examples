# function to generate a random combination of dispersal, filtering and competition parameter settings (uninformed prior assumed)
getRandomVals <- function(max_val)
{
  x <- runif(3, min = 0, max = 1);  
  x <- x / sum(x); #normalize to 1
  x <- x * max_val; #translate to integers
  x2 <- floor(x); #round to integers
  while(sum(x2) != max_val)
  {
  		a <- sample(1:3,size=1,prob=x-floor(x));
  		x2[a] <- x2[a]+1;
  }
  
  return(x2)
}


# function to randomly draw a particle depending on it's weight
getFromPrevious <- function(inds, ws, disps, filts, comps)
{
  index <- sample(x = inds, size = 1, replace = TRUE, prob = ws);
   
  output <- c(disps[index], filts[index], comps[index]);
  return(output); 
}

# function to calculate the weight of a particle
calculateWeight <- function(params, target, sigM, dispVals, filtVals, compVals, numParticles, W)
{
  sum = 0.0;
  vals <- c();
  numP <- length(dispVals);
      
  for(i in 1:numP)
  { 
    prevP <- c(dispVals[i], filtVals[i], compVals[i]);
    diff <- params[target] - prevP[target];
    vals[i] <- W[i] * dnorm(diff, mean = 0, sd = sigM);   
  }
  return(1 / sum(vals));
}

# normalize all the weights of all particles such that they sum to 1
normalizeWeights <- function(x)
{
  sum_x <- sum(x)    
  x <- x / sum_x;
  return(x);
}

# function to randomly change the contribution of one of the processes:
perturb <- function(p, sigM)
{
  params <- p;
  max_number <- sum(p);
  numbers <- 1:3
  
  x <- sample(numbers, 3, replace=F);
  
  oldVal <- params[x[1]];
  
  params[x[1]] <- round(params [x[1]] + rnorm(1, mean = 0, sd = sigM), 0);
  params[x[1]] <- max(0, params[x[1]]);
  params[x[1]] <- min(max_number, params[x[1]]);  

  diff <- params [x[1]] - oldVal;  
  
  params[x[2]] <- params[x[2]] - diff;
  params[x[2]] <- max(0, params[x[2]]);
  params[x[2]] <- min(max_number, params[x[2]]);

  params[x[3]] <- max_number - (params[x[1]] + params[x[2]]);
  params[x[3]] <- max(0, params[x[3]]);
  params[x[3]] <- min(max_number, params[x[3]]);

  return(c(params, x[1]));
}

# function to calculate the fit
calculateDistance <- function(rich, even, div, opt_diff, obs, sd_vals)
{
  Fit_Rich <- (abs((rich - obs[,1])) / sd_vals[1])^2;
  Fit_Even <- (abs((even - obs[,2])) / sd_vals[2])^2;
  Fit_Div  <- (abs((div -  obs[,3])) / sd_vals[3])^2;
  Fit_optima <- (opt_diff / sd_vals[4])^2;
    
  Full_Fit <- Fit_Rich + Fit_Even + Fit_Div + Fit_optima;
    
  return(Full_Fit)
}


ABC_SMC <- function(numParticles, species_fallout, taxa, esppres, n_traits,
sd_vals, summary_stats, community_number, species, abundances, frequencies, stopRate, Ord)
{
  optimum <- summary_stats[, 4:(3 + n_traits)];

  dispVals <- 1:numParticles;
  filtVals <- 1:numParticles;
  compVals <- 1:numParticles;

  fits <- 1:numParticles;
  RichVec <- 1:numParticles;
  EveVec <- 1:numParticles;
  DivVec <- 1:numParticles;
  OptVec <- 1:numParticles;

  nextDisp <- dispVals;
  nextFilt <- filtVals;
  nextComp <- compVals;

  weights <- rep(1, numParticles);
  nextWeights <- rep(1, numParticles);
  indices <- 1:numParticles;

  sigma <- 1;
  t <- 1;

  f <- list.files(pattern = "particles_t=");
   if(length(f) > 0)
  {
    f <- mixedsort(f)    
    t1 <- 1+length(f);
    d <- read.table(f[length(f)], header = F);
    
   
    dispVals <- d[,1];
    filtVals <- d[,2];
    compVals <- d[,3];
    fits <-     d[,8];
    weights <-  d[,9];	

    t <- t1;
	if(d[numParticles,1] == numParticles)
    {
       d <- read.table(f[length(f)-1],header=F);
	   output <- list( DA = d[,1],HF = d[,2], LS = d[,3]);   
	   cat("Found previously finished run, loaded results from that run\n"); flush.console();
	   return(output);
    }
  } 

  # continuously sampling
  while(t < 50)  
  {
    cat("\nGenerating Particles for iteration\t", t, "\n")
    
    cat("0--------25--------50--------75--------100\n")
    
    cat("*"); flush.console(); 
    PRINT_FREQ <- 20; 

    numberAccepted <- 0;
    if(t != 1) weights <- normalizeWeights(weights);
    threshold <- 200 * exp(-0.5 * t)
      
    stop_iteration <- 0;
    changed <- 1;

    tried <- 1;

    while(numberAccepted < numParticles)
    {
      params <- c(species_fallout, 0, 0);
      
      # get a parameter combination   
      if(t == 1)
      {
        params <- getRandomVals(species_fallout);
      }else{
        params <- getFromPrevious(indices, weights, dispVals, filtVals, compVals);
        params <- perturb(params, sigma);

        # we need to know which parameter was perturbed, 
        # to be able to calculate its weight later
        changed <- params[4]; 
        params <- params[1:3];
            
        if(sum(params) > species_fallout)
        { 
          print("too much params after perturb!"); flush.console(); break;
        } 
      }
      
      # total number of species in species pool 
      taxa <- length(abundances[1, ]) 
      allcommunities <- matrix(nrow = taxa)
      allcommunities <- STEPCAM(params, species, abundances, taxa,
      esppres, community_number, n_traits, species_fallout);

      # make traits (total species richness x number of traits) 
      # and community ((66 * permutations * communities) x total species richness) matrices
      traits <- as.data.frame(species[, c(2:(n_traits + 1))], row.names = c(1:taxa));
      communities <- as.data.frame(t(allcommunities), names = (1:taxa))
      names(communities) <- c(1:taxa)
      present_species <- as.vector(which(colSums(communities)>0))

      # calculate several measures of FD of modeled community
     # FD_output <- dbFD(traits[present_species, ], communities[, present_species], stand.x = F,messages=FALSE)
      FD_output <- strippedDbFd(Ord, communities[,present_species]);       
      FRic <- FD_output$FRic # FRic = functional richness (Villeger et al, 2008, Ecology)
      FEve <- FD_output$FEve # FEve = functional evenness (Villeger et al, 2008, Ecology)
      FDiv <- FD_output$FDiv # FDiv = functional diversity (Villeger et al, 2008, Ecology)
        
      trait_means <- c()
      for(i in 1:n_traits){
        # trait means of simulated community
        trait_means[i] <- mean(traits[present_species, i]) 
      }
      optimum_plus_trait_means <- rbind(optimum, trait_means);
      
      # calculate distance of trait mean simulated community from that of observed  
      mean_optimum <- dist(optimum_plus_trait_means);
      
      # (inverse) fit of model: euclidian distance of FD and trait mean values of
      # observed community from that of simulated 
      fit <- calculateDistance(FRic[[1]], FEve[[1]], FDiv[[1]], mean_optimum[1], summary_stats, sd_vals); 

      # function to accept / reject models based on the fit
      if(fit < threshold)
      {
        numberAccepted <- numberAccepted + 1
        nextDisp[numberAccepted] = params[1];
        nextFilt[numberAccepted] = params[2];
        nextComp[numberAccepted] = params[3];

        fits[numberAccepted] = fit;
        RichVec[numberAccepted] = FRic[[1]];
        EveVec[numberAccepted] = FEve[[1]];
        DivVec[numberAccepted] = FDiv[[1]];
        OptVec[numberAccepted] = mean_optimum;

        if(t == 1) {nextWeights[numberAccepted] = 1}
        else {
          nextWeights[numberAccepted] = calculateWeight(params, changed, sigma,
          dispVals, filtVals, compVals, numParticles, weights);
        }

	      if ((numberAccepted)%%(numParticles / PRINT_FREQ) == 0) 
        {
          cat("**") ; flush.console(); 
        } 
      }
        
      tried <- tried + 1;
      if(tried > (1/stopRate) && tried > 5)
      {
        # do not check every particle if the acceptance rate is OK
        if(numberAccepted / tried < stopRate) {stop_iteration <- 1; break;}
      }              
    }
     
    # replace values
    dispVals <- nextDisp;
    filtVals <- nextFilt;
    compVals <- nextComp;
    weights <- nextWeights; 

    output <- cbind(dispVals, filtVals, compVals, RichVec, EveVec, DivVec, OptVec, fits, weights);
    file_name <- paste("particles_t=", t, ".txt", sep="", collapse = NULL);
    write.table(output, file_name, row.names = F, col.names = F);
     
    cat(" ", mean(dispVals), mean(filtVals), mean(compVals), "\t", "accept rate = ",
    numberAccepted / (tried-1), numberAccepted, tried, "\n");
    
    # and reset
    nextWeights <- rep(1,numParticles);
    nextDisp <- 1:numParticles;
    nextFilt <- 1:numParticles;
    nextComp <- 1:numParticles; 
    t <- t + 1;
    if(stop_iteration == 1){break;}
  }
  if(t >= 2) {
	d <- read.table(paste("particles_t=", t - 2, ".txt", sep="", collapse = NULL), header = F);
  } else {
    if(t >= 1) {
		d <- read.table(paste("particles_t=", t - 1, ".txt", sep="", collapse = NULL), header = F);
	} else {
		d <- read.table(paste("particles_t=", t, ".txt", sep="", collapse = NULL), header = F);
	}
  }
  output <- list( DA = d[, 1], HF = d[, 2], LS = d[, 3]);
  return(output);
}


