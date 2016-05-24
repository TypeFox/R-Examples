#### 1- Computation of the spatial model ####

calcPotts <- function(W_SR, sample, rho, prior = TRUE, site_order = NULL, 
                      W_LR = NULL, nbGroup_min = 100, coords = NULL, distance.ref = NULL, threshold = 0.1, multiV = TRUE, 
                      iter_max = 200, cv.criterion = 0.005, verbose = 2){
  
  neutre <- 0 
  if(is.data.frame(coords)){ coords <- as.matrix(coords) }
  if(is.data.frame(sample)){ sample <- as.matrix(sample) }
  
  # initPackage(package = "spam", method = "calcPotts")
  # initPackage(package = "Matrix", method = "calcPotts")
  validInteger(iter_max,  validLength = 1, min = 1, method = "calcPotts")
  
  #### preparation
  n <- nrow(as.matrix(sample))
  p <- ncol(sample)  
  
  validW(W_SR, name = "W_SR", row.norm = TRUE , method = "calcPotts")
  validDim_matrix(W_SR, value2 = c(n, n), name1 = "W_SR", name2 = NULL, method = "calcPotts")
  
  if(!is.null(W_LR)){
    
    test.regional <- TRUE     
    validDim_vector(value1 = rho, value2 = list(2), name1 = "rho", name2 = NULL, type = "length", method = "calcPotts") 
    validDim_matrix(W_LR, value2 = c(n, n), name1 = "W_LR", name2 = NULL, method = "calcPotts")
    
    if(is.null(distance.ref)){
      distance.ref <- unlist(lapply(W_LR, function(x){sort(unique(x@x))}))
      distance.ref <- unique(distance.ref)
      distance.ref <- distance.ref[distance.ref != 0]
    }
    
    if(any(W_LR@x %% 1 > 0) || any(W_LR@x >= length(distance.ref)) ){
      W_LR@x <- findInterval(x = W_LR@x, vec = distance.ref) - 1
    }
    
    if(is.null(coords)){
      stop("calcPotts : argument \'coords\' cannot be NULL when regional regularization is ask \n")
    }

    if(is.data.frame(coords)){
      coords <- as.matrix(coords)
    }
    
    validDim_matrix(value1 = coords, value2 = n, name1 = "coords", name2 = n, type = "nrow", method = "calcPotts")
    
  }else{
    test.regional <- FALSE
    distance.ref <- 1:10 # evite un bug de la fonction qui veut un NumericVector et non NULL
    W_LR <- Matrix::Matrix(NA)
    coords <- matrix(NA, 1, 1)
  }
  
  if(test.regional == TRUE && rho[2] == 0){
    test.regional <- FALSE
  }
  
  
  #### initialisation  
  Wpred <- matrix(0, nrow = n, ncol = p)
  if(is.null(site_order)){
    if(verbose > 1){cat("NOTE : specifying argument \'site_order\' would speed up the execution of the function \n")}
    site_order <- -1
  }
  
  #### estimation
  
  res <- calcPotts_cpp(W_SR = spam::t(W_SR), W_LR = W_LR, sample = sample, 
                       rho = rho, coords = coords, site_order = site_order, iter_max = iter_max, cv_criterion = cv.criterion, 
                       test_regional = test.regional, distance_ref = distance.ref, threshold = threshold, neutre = neutre, 
                       nbGroup_min = nbGroup_min, multiV = multiV, last_vs_others = FALSE, prior = prior, type_reg = FALSE, 
                       verbose = verbose)
  
  res$predicted <- matrix(unlist(res$predicted), nrow = n, ncol = p, byrow = FALSE)
  res$radius <- list()
  res$Vregional <- matrix(NA, ncol = p, nrow = n)
  
  #### export
  if(verbose > 0 && res$cv == FALSE){
    cat("WARNING : the algorithm has not converged  \n", 
        "change argument \'iter_max\' to a higher value \n")
  }
  if(test.regional == TRUE && verbose > 1){
    
    for(iter_p in 1:p){
      
      resV <- calcMultiPotential_cpp(W_SR = W_SR, W_LR = W_LR, sample = res$predicted[,iter_p], 
                                     threshold = threshold, coords = coords, distance_ref = distance.ref, 
                                     nbGroup_min = nbGroup_min, multiV = multiV, neutre = neutre)
      res$radius[[iter_p]] <- resV$Radius
      res$Vregional[,iter_p] <- resV$V
      
      cat(" Group ", iter_p, " - radius : ", paste(round(res$radius[[iter_p]][res$radius[[iter_p]] != -1], optionsMRIaggr("digit.result")), collapse = " "))
      nb_m1 <- sum(res$radius[[iter_p]] == -1)
      if(nb_m1 > 0){cat(" -1 (#", nb_m1, ")")}
      cat("\n")
    }
    
  }
  
  return(res)
}

calcMultiPotential <- function(W_SR, W_LR, distance.ref, sample, coords, 
                               threshold = 0.01, nbGroup_min = 100, check.args = TRUE, verbose = TRUE){
  neutre <- 0 
  multiV <- TRUE
  
  # tests preliminaires
  if(check.args == TRUE){
    
    validW(W_SR, name = "W_SR", row.norm = TRUE, method = "calcMultiPotential")
    validW(W_LR, name = "W_LR", row.norm = FALSE, method = "calcMultiPotential")
    validInteger(nbGroup_min, validLength = 1, min = 0, method = "calcMultiPotential")
    validNumeric(threshold, validLength = 1, min = 0, method = "calcMultiPotential")    
    validClass(value = sample, validClass = "vector", superClasses = TRUE, method = "calcMultiPotential")
    
    n <- length(sample)
    
    if(is.data.frame(coords)){coords <- as.matrix(coords)}
    
    validClass(value = coords, name = "coords", validClass = "matrix", superClasses = TRUE, method = "calcMultiPotential")
    
    validDim_matrix(value1 = W_SR, value2 = c(n,n), name1 = "W_SR", name2 = NULL, type = "both", method = "calcMultiPotential")
    validDim_matrix(value1 = W_LR, value2 = c(n,n), name1 = "W_SR", name2 = NULL, type = "both", method = "calcMultiPotential")
    validDim_matrix(value1 = coords, value2 = n, name1 = "coords", name2 = NULL, type = "nrow", method = "calcMultiPotential")
    
    if(length(distance.ref) < 2 || any(distance.ref <= 0)){
      stop("calcMultiPotential : wrong specification of \'distance.ref\' \n", 
           "\'distance.ref\' must contains positive elements and be at least of length 2 \n", 
           "proposed \'distance.ref\' : ", paste(distance.ref, collapse = " "), "\n")
    }
    
    if(any(W_LR@x %% 1 > 0) || any(W_LR@x >= length(distance.ref)) ){
      
      if(verbose){cat("segment x slot of argument \'W_LR\' according to argument \'distance.ref\' \n")}
      W_LR@x <- findInterval(x = W_LR@x, vec = distance.ref) - 1
      
    }
    
    if(any(W_LR@x < 0)){
      stop("calcMultiPotential : wrong specification of \'W_LR\' \n", 
           "W_LR@x must contain integers between 0 and ", length(distance.ref) - 1, " \n")
    }
    
  }
  
  #### calcul du potentiel
  res <- calcMultiPotential_cpp(W_SR = W_SR, W_LR = W_LR, sample = sample, 
                                threshold = threshold, coords = coords, distance_ref = distance.ref, nbGroup_min = nbGroup_min, 
                                multiV = multiV, neutre = neutre)
  
  
  if(verbose == TRUE){
    cat("radius : ", paste(round(res$Radius, optionsMRIaggr("digit.result")), collaspse = " "), "\n")
  }
  
  return(res$V)
}


#### 2- estimation of the spatial parameter ####
calcPottsParameter <- function(Y, W_SR, coords = NULL, range = NULL, method = "MF", verbose = 3, ...){
  
  args <- list(...)
  names.args <- names(args)
  
  
  #### initialisation
  validCharacter(method, validLength = 1, validValues = c("MF", "Lvfree"), refuse.NULL = TRUE, method = "calcPottsParameter")
  
  if(missing(W_SR) || is.null(W_SR)){
    if(verbose > 0){cat("Computation of the local neighborhood matrix \n")}
    validClass(coords, validClass = "data.frame", superClasses = FALSE, method = "calcPottsParameter")
    
    resW <- calcW(coords, range = range[1], upper = NULL, row.norm = TRUE, calcBlockW = (method == "Lvfree"))
    W_SR <- resW$W
    if(method == "Lvfree"){
      site_order <- unlist(resW$blocks$ls_groups) - 1
    }
  }
  
  if(is.data.frame(Y)){
    Y <- as.matrix(Y)
  }
  
  validW(W_SR, name = "W_SR", row.norm = TRUE, method = "calcPottsParameter")
  
  #### Mean Field estimation
  if(method == "MF"){
    if(verbose > 0){cat("Estimation using mean field approximation \n")}
    
    #### tests
    validCharacter(value = names.args, name = "...", validLength = NULL, 
                   validValues = c("rho_max", "prior_prevalence", "W_LR", "distance.ref", "threshold", "nbGroup_min", "coords", "regionalGroups", "multiV"), 
                   refuse.NULL = FALSE, method = "calcPottsParameter")
    
    if("W_LR" %in% names.args && !is.null(args$W_LR)){
      test.regional <- TRUE
      if("regionalGroups" %in% names.args){
        validCharacter(args$regionalGroups, "regionalGroups", validLength = 1, validValues = c("last_vs_others", "each"), refuse.NULL = TRUE, method = "rhoMF")
      }
    }else{
      test.regional <- FALSE
    }
    
    if(test.regional == TRUE){
      
      validW(args$W_LR, name = "W_LR", row.norm = FALSE, method = "rhoMF")
      
      if("distance.ref" %in% names.args == FALSE || is.null(args$distance.ref)){
        args$distance.ref <- unlist(lapply(args$W_LR, function(x){sort(unique(x@x))}))
        args$distance.ref <- unique(args$distance.ref)
        args$distance.ref <- args$distance.ref[args$distance.ref != 0]
      }
      
      if(any(args$W_LR@x %% 1 > 0) || any(args$W_LR@x >= length(args$distance.ref)) ){
        args$W_LR@x <- findInterval(x = args$W_LR@x, vec = args$distance.ref) - 1
      }
      
      if("coords" %in% names.args == FALSE || is.null(args$coords)){
        stop("calcPottsParameter : argument \'coords\' must not be NULL if argument \'W_LR\' is specified \n")
      }
      
      if(is.data.frame(args$coords)){
        args$coords <- as.matrix(args$coords) 
      }
      
    }
    
    #### estimation 
    rho <- rhoMF(Y = Y, W_SR = W_SR, test.regional = test.regional, ...)
    
    if(verbose > 0){
      cat("Estimated rho : ", paste(round(rho, optionsMRIaggr("digit.result")), collpase = " "), "\n", sep = "")
    }
  }
  
  #### Likelihood Free estimation
  if(method == "Lvfree"){
    if(verbose > 0){cat("Estimation using a likelihood free method \n")}
    
    #### tests
    validCharacter(value = names.args, name = "...", validLength = NULL, 
                   validValues = c("rho_max", "rho_init", "sd_rho_init", "site_order", "dprior", 
                                   "epsilon", "update_epsilon", "n.sample", "iter_max", "burn_in", 
                                   "thin", "n.chain", "verbose", "export.coda", "cpus"), 
                   refuse.NULL = TRUE, method = "calcPottsParameter")
    
    if("n.sample" %in% names.args && args$n.sample < 50){
      stop("calcPottsParameter : \'n.sample\' too low \n", 
           "\'n.sample\' must be al least 50 \n", 
           "proposed \'n.sample\' = ", args$n.sample, "\n")    
    }
    
    
    if(length(args$site_order) == 1 && args$site_order == TRUE){
      args$site_order <- unlist(calcBlockW(W_SR, verbose = verbose)$ls_groups) - 1
    }
    
	n <- nrow(Y)
    validInteger(value = args$site_order, name = "args$site_order", validLength = n, 
                 min = 0, max = n - 1, refuse.NA = TRUE, refuse.NULL = FALSE, refuse.duplicates = TRUE, method = "calcPottsParameter")
    
    
    if(is.null(site_order)){	
      if(verbose == TRUE){cat("NOTE : specifying argument \'site_order\' would speed up the execution of the function \n")}
      args$site_order <- -1	 
    }
    
    #### estimation
    
    rho <- rhoLvfree(Y = Y, W_SR = W_SR, verbose = verbose, ...)
    
    if(verbose > 0 && sum(!is.na(rho)) > 0){
      cat("Estimated rho (by each chain) : ", paste(unlist(lapply(rho, 
                                                                  function(x){paste(round(stats::median(x), digits = optionsMRIaggr("digit.result")), 
                                                                                    " [", round(stats::quantile(x, 0.025), digits = optionsMRIaggr("digit.result")), 
                                                                                    " ; ", round(stats::quantile(x, 0.975), digits = optionsMRIaggr("digit.result")), "]",
                                                                                    sep = "")})
      ), collpase = " "), "\n", sep = "")
    }
    
  }
  
  #### export
  return(rho)
}

rhoMF <- function(Y, W_SR, rho_max = 50, prior_prevalence = TRUE, 
                  test.regional = FALSE, W_LR, distance.ref, coords, threshold = 0.1, nbGroup_min = 100, 
                  regionalGroups = "last_vs_others", multiV = TRUE){
  
  neutre <- 0.5
  
  #### preparation
  n <- nrow(as.matrix(Y))
  p <- ncol(Y)
  prevalence <- colMeans(Y)
  
  #### estimation du rho
  
  # calcul du potentiel local
  U <- matrix(NA, nrow = n, ncol = p)
  for(iter_p in 1:p){
    U[,iter_p] <- as.vector( get("%*%", envir = asNamespace("spam"))(W_SR, Y[,iter_p]) )
  }
  
  # calcul du potentiel regional
  V <- matrix(NA, n, p)
  
  if(test.regional == TRUE){
    
    if(regionalGroups == "last_vs_others"){
      V[,p] <- calcMultiPotential_cpp(W_SR = W_SR, W_LR = W_LR, sample = Y[,p], threshold = threshold, 
                                      coords = coords, distance_ref = distance.ref, nbGroup_min = nbGroup_min, 
                                      multiV = multiV, neutre = neutre)$V
      
      V[,-p] <- matrix(1 - V[,p], n, p - 1, byrow = FALSE)
      
    }else if(regionalGroups == "each"){
      
      for(iter_regionalGroup in 1:p){
        V[,iter_regionalGroup] <- calcMultiPotential_cpp(W_SR = W_SR, W_LR = W_LR, sample = Y[,iter_regionalGroup], threshold = threshold, 
                                                         coords = coords, distance_ref = distance.ref, nbGroup_min = nbGroup_min, 
                                                         multiV = multiV, neutre = neutre)$V
      }
      
    }
    
  } 
  
  ### fonction objectif
  E <- matrix(NA, n, p)
  PML <- function(rho){
    ###  estimation conjointe
    # Y : valeurs observees (binaire, eventuellement entre 0 et 1)
    # U = WY : moyenne des valeur des voisins proches
    # E : matrice de NA n, p
    # poids_prevalence : vecteur de poids des observations
    
    E <-   rho * U
    K <- rowSums(exp(E))
    lv <- sum( rowSums( (E-log(K)) * Y) * poids_prevalence )
    return(-lv)
    
  }
  dPML <- function(rho){
    ###  estimation conjointe
    # Y : valeurs observees (binaire, eventuellement entre 0 et 1)
    # U = WY : moyenne des valeur des voisins proches
    # V : densite de voisins regionale
    # E : matrice de NA n, p
    # poids_prevalence : vecteur de poids des observations
    
    E <-   rho * U
    K <- rowSums(U * exp(E)) / rowSums(exp(E))
    dPML <- apply(U, 2, function(x){x - K})
    dPML <- sum( rowSums(dPML * Y) * poids_prevalence )
    return(-dPML)
  }
  
  PML_regional <- function(rho){
    ###  estimation conjointe
    # Y : valeurs observees (binaire, eventuellement entre 0 et 1)
    # U = WY : moyenne des valeur des voisins proches
    # V : densite de voisins regionale
    # E : matrice de NA n, p
    # poids_prevalence : vecteur de poids des observations
    
    E <- rho[1] * U + rho[2] * V
    
    K <- rowSums(exp(E))
    lv <- sum( rowSums( (E-log(K)) * Y) * poids_prevalence )
    return(-lv)
    
  }
  
  dPML_regional <- function(rho){
    ###  estimation conjointe
    # Y : valeurs observees (binaire, eventuellement entre 0 et 1)
    # U = WY : moyenne des valeur des voisins proches
    # V : densite de voisins regionale
    # E : matrice de NA n, p
    # poids_prevalence : vecteur de poids des observations
    
    E <- rho[1] * U + rho[2] * V
    
    K1 <- rowSums(U * exp(E)) / rowSums(exp(E))
    dPML1 <- apply(U, 2, function(x){x - K1})
    dPML1 <- sum( rowSums(dPML1 * Y) * poids_prevalence )
    
    K2 <- rowSums(V * exp(E)) / rowSums(exp(E))
    dPML2 <- apply(V, 2, function(x){x - K2})
    dPML2 <- sum( rowSums(dPML2 * Y) * poids_prevalence )
    
    return(-c(dPML1, dPML2))
  }
  
  
  if(prior_prevalence == FALSE){
    poids_prevalence <- as.matrix(Y) %*% cbind(1 / prevalence)
  }else{
    poids_prevalence <- 1
  }
  
  if(test.regional == TRUE){
    rho <- optim(par = c(0, 0), method = "L-BFGS-B", lower = c(0, 0), upper = c(rho_max, rho_max), 
                 fn = PML_regional, gr = dPML_regional)$par
  }else{
    rho <- optim(par = 0, method = "L-BFGS-B", lower = 0, upper = rho_max, 
                 fn = PML, gr = dPML)$par
  }
  
  #### export ####
  return(rho)
}


rhoLvfree <- function(Y, W_SR, rho_max = 10, rho_init = "PML", sd_rho_init = "PML", site_order = NULL, # to be initialized
                      dprior = function(x){stats::dunif(x, 0, rho_max)}, epsilon = 0.001, update_epsilon = 1,                      
                      n.sample = 2500, iter_max = 3, burn_in = round(0.25 * n.sample), thin = 1, n.chain = 1, 
                      verbose = 3, cpus = 1, export.coda = TRUE){
  
  n <- nrow(Y)
  G <- ncol(Y)
  
  ## initialisation
  if(length(rho_init == 1) && rho_init == "PML"){rho_init <- rhoMF(Y = Y, W_SR = W_SR, rho_max = rho_max)}    
  if(length(rho_init) == 1){rho_init <- rep(rho_init, n.chain)}
  if(length(sd_rho_init == 1) && sd_rho_init == "PML"){
    sd_rho_init <- if(n < 500){0.1 * rho_init}else{if(n < 5000){0.05 * rho_init}else{0.025 * rho_init}}
    sd_rho_init[sd_rho_init < 0.25] <- 0.25
  }    
  if(length(sd_rho_init) == 1){sd_rho_init <- rep(sd_rho_init, n.chain)}
  Uobs <- sum(Y * get("%*%", envir = asNamespace("spam"))(W_SR, Y) )
  tW_SR <- spam::t(W_SR)
  
  ls.MCMC_rho <- list()
  
  ## fct
  warper_chain <- function(iter_chaine, rho, sd_rho, burn_in){
    
    MCMC_rho <- rep(NA, n.sample)
    epsilon_chain <- epsilon
    burn_in_chain <- burn_in
    acceptation <- rep(NA, n.sample)
    test.burn_in <- FALSE
    
    if(verbose > 0){cat("> chain ", iter_chaine, " (rho init = ", rho, " ; sd init = ", sd_rho, " ) \n", sep = "")}
    
    ## echantillonnage 
    for(iter in 1:n.sample){
      # initialisation
      rho_test <- rtnorm(1, mean = rho, sd = sd_rho, lower = 0, upper = rho_max)
      
      # simulation du champ
      gibbs <- simulPottsFast_cpp(W_i = W_SR@i, W_p = W_SR@p, W_x = W_SR@x, 
                                  site_order = site_order, sample = Y, rho = rho_test, n = n, p = G, 
                                  iter_nb = iter_max)       
      
      # evaluation de la statistique exhaustive
      U <- sum(gibbs * get("%*%", envir = asNamespace("spam"))(W_SR, gibbs) )
      
      # maj du beta
      if( abs(U - Uobs) < epsilon_chain * Uobs){
        
        ratio <- dprior(rho_test) * dtnorm(rho, mean = rho_test, sd = sd_rho, lower = 0, upper = rho_max) / (dprior(rho) * dtnorm(rho_test, mean = rho, sd = sd_rho, lower = 0, upper = rho_max))
        if(ratio > stats::runif(1, 0, 1)){
          rho <- rho_test
          acceptation[iter] <- rho_test         
        }
      }
      
      # gestion de l historique et parametres
      MCMC_rho[iter] <- rho
      
      # affichage
      if(iter %% round(0.025 * n.sample) == 0){
        if(verbose > 1){ cat("*") }
        
        if(test.burn_in == TRUE){
          sd_rho <- sd(acceptation, na.rm = TRUE) 
          
          if(update_epsilon > 1){ # reduce epsilon
            epsilon_chain <- epsilon_chain / update_epsilon
          }else if(update_epsilon < 1){ # set acceptance rate : update_epsilon should correspond to the target acceptance rate
            epsilon_chain <- epsilon_chain * (1 + (update_epsilon - length(na.omit(acceptation)) / iter))
          }
          
          if(verbose == 3){ cat(" sd rho = ", round(sd_rho, optionsMRIaggr("digit.result")), 
                                " - acceptance = ", round(length(na.omit(acceptation)) / iter, optionsMRIaggr("digit.percentage")), 
                                " - epsilon = ", round(epsilon_chain, optionsMRIaggr("digit.epsilon")), "\n") }
        }
      }
      
      if(iter == burn_in_chain){
        
        sd_rho_test <- sd(acceptation, na.rm = TRUE)
        
        if(is.na(sd_rho_test) || sd_rho_test == 0){
          if(verbose > 1){
            warning("rhoLvfree : \'burn_in_chain\' too short \n", 
                    "no rho accepted, null variance \n")
          }
          burn_in_chain <- max(burn_in_chain * 2, n.sample)
        }else{
          sd_rho <- sd_rho_test
          test.burn_in <- TRUE
          if(verbose == 3){
            cat(" sd rho = ", round(sd_rho, optionsMRIaggr("digit.result")), 
                " - acceptance = ", round(length(na.omit(acceptation)) / iter, optionsMRIaggr("digit.percentage")), 
                " - epsilon = ", round(epsilon_chain, optionsMRIaggr("digit.epsilon")), "\n") }
        
        }
      }
    }
    
    # burn in
    if(test.burn_in && burn_in_chain < n.sample){
      MCMC_rho <- MCMC_rho[seq(from = burn_in_chain + 1, to = n.sample, by = thin)]
    }else{
      MCMC_rho <- NA
    }
    if(verbose > 1){cat("\n")}
    
    return(MCMC_rho)
  }
  
  ## loop
  if(cpus == 1){
    ls.MCMC_rho <- lapply(1:n.chain, 
                          function(x){warper_chain(x, rho = rho_init[x], sd_rho = sd_rho_init[x], burn_in = burn_in)}
    )
  }else{
    initPackage(package = "package", argument = "\'cpus\' > 1", tryAttach = TRUE, method = "snowfall")
    initPackage(package = "snowfall", argument = "\'cpus\' > 1", method = "snowfall")
    
    snowfall::sfInit(parallel = TRUE, cpus = cpus)
    snowfall::sfLibrary( "MRIaggr" )
    verbose <- 0
    snowfall::sfExport("burn_in", "dprior", "epsilon", "G", "iter_max", 
                       "n", "n.sample", "rho_init", "rho_max", 
                       "sd_rho_init", "site_order", "thin", "tW_SR", "Uobs", 
                       "update_epsilon", "W_SR", "warper_chain", "Y")
    
    ls.MCMC_rho <- snowfall::sfClusterApplyLB(1:n.chain, 
                                              function(x){warper_chain(x, rho = rho_init[x], sd_rho = sd_rho_init[x])}
    )
    
    snowfall::sfStop()
  }
  
  # export en coda
  if(export.coda == TRUE){
    initPackage(package = "coda", argument = "\'export.coda\' = TRUE", method = "rhoLvfree")
    size_tempo <- min(unlist(lapply(ls.MCMC_rho, length)))
    ls.MCMC_rho <- lapply(ls.MCMC_rho, function(x){tail(x, size_tempo)})
    ls.MCMC_rho <- lapply(ls.MCMC_rho, coda::as.mcmc)
    
    ls.MCMC_rho <- coda::as.mcmc.list(ls.MCMC_rho)
  }
  
  return(ls.MCMC_rho)
  
}


#### 3- simulation of a spatial field ####

simulPotts <- function(W, G, rho, initialization = NULL, site_order = NULL, iter_max = 2500, 
                       cv.criterion = 0, fast = FALSE, verbose = TRUE){
  
  #### test
  
  # initialization
  if(!is.null(initialization)){
    
    validClass(value = initialization, name = "initialization", validClass = c("NULL","matrix"), superClasses = FALSE, method = "simulPotts")
    
    G <- ncol(initialization)
    if(G < 2){
      stop("simulPotts : incorrect specification of \'initialization\' \n", 
           "\'initialization\' must have at least two columns \n", 
           "ncol(initialization) : ", G, "\n")
    }
    
    n <- nrow(initialization)
    validDim_matrix(value1 = W, value2 = n, name1 = "W", name2 = "initialization", type = "nrow", method = "simulPotts")
    
  }else{
    validInteger(G, validLength = 1, min = 2, method = "simulPotts")
    n <- nrow(W)
    initialization <- t(stats::rmultinom(n, size = 1, prob = rep(1 / G, G)))
  }
  
  # site order
  if(length(site_order) == 1 && site_order == TRUE){
    site_order <- unlist(calcBlockW(W, verbose = verbose)$ls_groups) - 1
  }
  
  validInteger(value = site_order, validLength = n, 
               min = 0, max = n - 1, refuse.NULL = FALSE, refuse.duplicates = TRUE, method = "calcPottsParameter")
  
  
  if(is.null(site_order)){	
    if(verbose == TRUE){cat("NOTE : specifying argument \'site_order\' would speed up the execution of the function \n")}
    site_order <- -1	 
  }
  
  # miscellaneous
  validInteger(iter_max, validLength = 1, min = 0, method = "simulPotts")
  validNumeric(cv.criterion, name = "cv.criterion", validLength = 1, min = 0, max = 1, method = "simulPotts")
  validLogical(verbose, validLength = 1, method = "simulPotts")
  validLogical(fast, validLength = 1, method = "simulPotts")
  
  if(verbose == TRUE && fast == TRUE){
    stop("simulPotts : \'fast\' cannot be set to TRUE is \'verbose\' is TRUE \n")
  }
  if(cv.criterion > 0 && fast == TRUE){
    stop("simulPotts : \'fast\' cannot be set to TRUE is \'cv.criterion\' is not 0 \n")
  }
  
  #### fct
  
  if(fast == TRUE){
    res <- simulPottsFast_cpp(W_i = W@i, W_p = W@p, W_x = W@x, site_order = site_order, 
                              sample = initialization, rho = rho, n = n, p = G, 
                              iter_nb = iter_max)
    
    return(list(simulation = res, 
                V = NA, 
                cv.criterion = NA, 
                cv = NA))
  }else{
    return(simulPotts_cpp(W_SR = W, W_LR = Matrix::sparseMatrix(i = 1, j = 1, x = 1), 
                          sample = initialization, coords = matrix(0), site_order = site_order,  
                          rho = c(rho, 0), distance_ref = c(0), 
                          iter_max = iter_max, cv_criterion = cv.criterion, regional = FALSE, 
                          verbose = verbose)
    )
    
  }
  
}    


#### 4- valid functions ####
validW <- function(value, name, row.norm, method){
  
  if("dgCMatrix" %in% is(value) == FALSE){
    stop("simulPotts : wrong specification of \'", name, "\' \n", 
         "\'", name, "\' must be a sparse matrix of class \'dgCMatrix\' \n", 
         "is(", name, ") : ", paste(is(value), collapse = " "), "\n", 
         "use the calcW function to compute ", name, " \n")
  }
  
  if(any(value@x < 0)){
    stop("simulPotts : wrong specification of \'", name, "\' \n", 
         "\'", name, "\' must only contain positive values \n", 
         "use the calcW function to compute ", name, " \n")
  }
  
  if(row.norm){
    WrowSum <- spam::rowSums(value)
    
    if(any(abs(WrowSum[WrowSum > 0] - 1) > 0.0001)){
      warning("simulPotts : \'", name, "\' is not row normalized \n", 
              "it may lead to an incorrect simulation \n", 
              "use the calcW function to compute ", name, " with argument row.norm = TRUE")
    }
  }
  
}
