#### 1- calc functions ####
# Otsu non retenu car si deux groupes 1 et 0
# Otsu prend une observation du groupe 1 (mean f^2  = 1) et le reste dans l'autre groupe (mean f^2 > 0)
# plutot que tout le groupe 1 (mean f^2  = 1) vs tout l autre groupe (mean f^2 = 0)

calcCriteriaGR <- function(contrast, groups, W = NULL, sigma = NULL, breaks = NULL, rm.warning = TRUE, 
                           criterion.transition = FALSE, criterion.sdfront = FALSE, criterion.entropy = TRUE, criterion.Kalinsky = TRUE, criterion.Laboure = TRUE){
  
  #### initialization
  # initPackage(package = "spam", method = "calcCriteriaGR")
  # initPackage(package = "Matrix", method = "calcCriteriaGR")
  
  n <- length(contrast) 
  
  criteria <- c("transition", "sdfront", "entropy", "Kalinsky", "Laboure")[c(criterion.transition, criterion.sdfront, criterion.entropy, criterion.Kalinsky, criterion.Laboure)]
  df.criterion <- data.frame(matrix(NA, ncol = 1 + length(criteria), nrow = 1))
  names(df.criterion) <- c("sigma", criteria)
  
  if(!is.null(sigma)){df.criterion$sigma <- sigma}
  
  indexT <- which(groups)
  n.T <- length(indexT)
  indexF <- which(groups == FALSE)
  n.F <- length(indexF)
  
  if(length(indexT) < 2 || length(indexF) < 2){
    if(rm.warning == TRUE){
      warning("calcCriteriaGR :  too small groups \n", 
              "sum(groups == TRUE) : ", sum(groups == TRUE), "\n", 
              "sum(groups == FALSE) : ", sum(groups == FALSE), "\n")
    }
    return(df.criterion)
  }
  
  #### criterion sur la frontiere ####
  if(criterion.transition == TRUE || criterion.sdfront == TRUE){
    
    # identification de la frontiere
    nb.voisT <- Matrix::rowSums(W[indexT, indexF] > 0)
    index.frontiereT <- indexT[which(nb.voisT > 0)]
    n.frontiere <- length(index.frontiereT)
    
    # iteration sur les points de  la frontiere
    if(criterion.transition == TRUE){
      df.criterion$transition <- 0 
    }
    if(criterion.sdfront == TRUE){
      df.criterion$sdfront <- 0
      nb_pointsF <- 0
    }
    
    for(iter_frontiere in 1:n.frontiere){
      indexiter_frontiere <- index.frontiereT[iter_frontiere]
      indexiter_voisF <- indexF[which(W[indexiter_frontiere, indexF] > 0) ]
      
      if(criterion.transition == TRUE){
        df.criterion$transition <- df.criterion$transition + sum(abs(contrast[indexiter_frontiere]-contrast[indexiter_voisF]))
        nb_pointsF <- nb_pointsF + length(indexiter_voisF)
      }
      
      if(criterion.sdfront == TRUE){
        df.criterion$sdfront <- df.criterion$sdfront + sum( ( contrast[indexiter_voisF]-contrast[indexiter_frontiere] ) )^2 
      }
    }
    
    if(criterion.transition == TRUE){df.criterion$transition <- df.criterion$transition / nb_pointsF}
    if(criterion.sdfront == TRUE){df.criterion$sdfront <- sqrt(df.criterion$sdfront / n.frontiere)}
  }
  
  #### criterion sur l'entropy ####
  if(criterion.entropy == TRUE){
    
    newbreaks <- graphics::hist(contrast, breaks = breaks, plot = FALSE)$breaks
    RF <- graphics::hist(contrast[indexF], breaks = newbreaks, plot = FALSE)$counts / length(indexF)
    RT <- graphics::hist(contrast[indexT], breaks = newbreaks, plot = FALSE)$counts / length(indexT)
    
    RF <- RF[RF > 0] # car x * log(x) tend vers 0 si x tend vers 0
    RT <- RT[RT > 0]
    
    df.criterion$entropy <- - sum(RF * log(RF)) - sum(RT * log(RT))
  } # http://ij-plugins.sourceforge.net/plugins/segmentation/Maximum_Entropy_Thresholding.pdf
  
  #### criterion sur la region ####
  if(criterion.Kalinsky == TRUE || criterion.Laboure == TRUE){
    moyT <- mean(contrast[indexT])
    moyF <- mean(contrast[indexF])
    WT <- sum((contrast[indexT]-moyT)^2)
    WF <- sum((contrast[indexF]-moyF)^2)
  }
  
  if(criterion.Kalinsky == TRUE){
    
    moy <- mean(contrast)
    
    B <- (moy - moyT)^2 + (moy - moyF)^2 # ddf(B) = 2 - 1 = 1
    W <- (WT + WF) / (n.T + n.F-1) # ddf(W) = n - 1
    df.criterion$Kalinsky <- B / W
  }
  
  if(criterion.Laboure == TRUE){
    df.criterion$Laboure <- 1 / sqrt(WT + WF)
  }
  
  return(df.criterion)
}



calcGR <- function(contrast, W, seed, sigma_max, range = c(-Inf, +Inf), range.seed = c(-Inf, +Inf), 
                   breaks = 100, rescale = FALSE, iter_max = 100, sd.robust = FALSE,                   
                   keep.lower = FALSE, keep.upper = FALSE, verbose = TRUE,                   
                   history.sigma = FALSE, history.step = FALSE, history.front = FALSE){
  
  #### initialization and tests
  res.initGR <- initGR(contrast = contrast, W = W, seed = seed, sigma_max = sigma_max, range = range, range.seed = range.seed, 
                       breaks = breaks, iter_max = iter_max, sd.robust = sd.robust, keep.lower = keep.lower, keep.upper = keep.upper,
                       history.sigma = history.sigma, history.step = history.step, history.front = history.front,
                       rescale = rescale, verbose = verbose, method = "calcGR")
  
  contrast <- res.initGR$contrast
  breaks <- res.initGR$breaks 
  step <- res.initGR$step 
  seed <- res.initGR$seed 
  
  #### computation : algorithm 
  if(verbose == TRUE){cat("GR (sigma_max = ", sigma_max, ") : ", sep = "")}
  
  res.GR <- GRalgo(contrast = contrast, W = W, seed = seed, sigma_max = sigma_max, range = range, keep.lower = keep.lower, keep.upper = keep.upper, 
                  breaks = breaks, step = step, operator = if(sd.robust){"stats::mad"}else{"stats::sd"}, iter_max = iter_max,                    
                  history.sigma = history.sigma, history.step = history.step, history.front = history.front)
  
  if(verbose == TRUE){cat(length(res.GR$GR), " observations in ", res.GR$iter, " iterations \n", sep = "")}
  
  if(res.GR$iter == 0 || (res.GR$iter == iter_max && res.GR$test.id == FALSE)){
    cv <- FALSE
    if(verbose == TRUE){
      warning("calcGR : the algorithm has not converged \n", 
              "maximum number of iterations reached (iter_max = ", iter_max, ")\n")
    }    
  }else if(res.GR$test.break == TRUE){
    warning("calcGR : the convergence may be unsatisfactory \n")
  }
  
  #### export
  res.GR$breaks <- breaks
  
  if(history.front && res.GR$iter > 1){
    front.split <- lapply(strsplit(res.GR$history.front, split = ".", fixed = TRUE), 
                          function(x){c(as.numeric(x), rep(NA, res.GR$iter-length(x)))}
    )
    res.GR$Mfront <- matrix(unlist(front.split), ncol = res.GR$iter, nrow = length(res.GR$GR), byrow = TRUE)
  }else{
    res.GR$Mfront <- NULL
  }
  
  res.GR$breaks <- breaks
  return(res.GR)
  
}

calcSigmaGR <- function(contrast, W, seed, sigma, 
                        criterion.transition = FALSE, criterion.sdfront = FALSE, criterion.entropy = TRUE, criterion.Kalinsky = TRUE, criterion.Laboure = TRUE, 
                        verbose = TRUE, ...){
  
  
  #### pre initialization
  n <- length(contrast)
  
  sigma <- sort(sigma)
  
  dots.arguments <- list(...)
  names_dots.arguments <- names(dots.arguments)
   
  if("range" %in% names(dots.arguments) == FALSE){dots.arguments$range <- c(-Inf, +Inf)}
  if("range.seed" %in% names(dots.arguments) == FALSE){dots.arguments$range.seed <- c(-Inf, +Inf)}
  if("breaks" %in% names(dots.arguments) == FALSE){dots.arguments$breaks <- 100}
  if("rescale" %in% names(dots.arguments) == FALSE){dots.arguments$rescale <- FALSE}
  if("iter_max" %in% names(dots.arguments) == FALSE){dots.arguments$iter_max <- 100}
  if("sd.robust" %in% names(dots.arguments) == FALSE){dots.arguments$sd.robust <- FALSE}
  if("keep.lower" %in% names(dots.arguments) == FALSE){dots.arguments$keep.lower <- FALSE}
  if("keep.upper" %in% names(dots.arguments) == FALSE){dots.arguments$keep.upper <- FALSE}
  if("rm.warning" %in% names(dots.arguments) == FALSE){dots.arguments$rm.warning <- FALSE}
  
  #### tests
  if (optionsMRIaggr("checkArguments")) {
    
    validNames(dots.arguments, name = "...", 
               validValues = c("range", "range.seed", "breaks", "rescale", "iter_max", "sd.robust", "keep.lower", "keep.upper", "rm.warning"), 
               method = "calcSigmaGR")
    
    validNumeric(value = sigma, validLength = NULL, min = 0, refuse.duplicates = TRUE, method = "calcSigmaGR")
    validLogical(value = dots.arguments$rm.warning, validLength = 1, method = "calcCriteriaGR")
    validLogical(value = criterion.transition, validLength = 1, method = "calcCriteriaGR")
    validLogical(value = criterion.sdfront, validLength = 1, method = "calcCriteriaGR")
    validLogical(value = criterion.entropy, validLength = 1, method = "calcCriteriaGR")
    validLogical(value = criterion.Kalinsky, validLength = 1, method = "calcCriteriaGR")
    validLogical(value = criterion.Laboure, validLength = 1, method = "calcCriteriaGR")
    if(criterion.transition == FALSE && criterion.sdfront == FALSE && criterion.entropy == FALSE && criterion.Kalinsky == FALSE && criterion.Laboure == FALSE){
      stop("calcSigmaGR : all criteria must not be simultaneously FALSE \n")
    }
    
  }
  
  #### initialization
  res.initGR <- initGR(contrast = contrast, W = W, seed = seed, sigma_max = sigma[1], range = dots.arguments$range, range.seed = dots.arguments$range.seed, 
                       breaks = dots.arguments$breaks, iter_max = dots.arguments$iter_max, sd.robust = dots.arguments$sd.robust, keep.lower = dots.arguments$keep.lower, keep.upper = dots.arguments$keep.upper,
                       history.sigma = TRUE, history.step = FALSE, history.front = FALSE,
                       rescale = dots.arguments$rescale, verbose = verbose, method = "calcSigmaGR")
  
  contrast <- res.initGR$contrast
  breaks <- res.initGR$breaks 
  step <- res.initGR$step 
  seed.sigma <- res.initGR$seed 

  n.sigma <- length(sigma)
  
  criteria <- c("transition", "sdfront", "entropy", "Kalinsky", "Laboure")[c(criterion.transition, criterion.sdfront, criterion.entropy, criterion.Kalinsky, criterion.Laboure)]
  n.criteria <- length(criteria)
  
  df.criterion <- data.frame(matrix(NA, ncol = 2 + n.criteria, nrow = n.sigma))
  names(df.criterion) <- c("sigma", "n.obs", criteria)
  df.criterion$sigma <- sigma
  
  list.GR <- sapply(1:n.criteria, function(x){c(list(), NA)})
  names(list.GR) <- criteria
  
  best <- data.frame(matrix(-Inf, nrow = 1, ncol = n.criteria))
  names(best) <- criteria
  
  #### computation : loop 
  if(verbose == TRUE){cat("loop over ", n.sigma, " sigma : ", sep = "")}
  
  for(iter_sigma in 1:n.sigma){
    if(verbose == TRUE){cat("*")}
    
    ## GR    
    res.GR <- GRalgo(contrast = contrast, W = W, seed = seed.sigma, sigma_max = sigma[iter_sigma], range = dots.arguments$range, 
                    breaks = breaks, step = step, operator = if(dots.arguments$sd.robust){"stats::mad"}else{"stats::sd"}, iter_max = dots.arguments$iter_max, 
                    keep.lower = dots.arguments$keep.lower, keep.upper = dots.arguments$keep.upper,                    
                    history.sigma = TRUE, history.step = FALSE, history.front = FALSE)
    
    ## compute criteria
    df.criterion$n.obs[iter_sigma] <- length(res.GR$GR)
    df.criterion[iter_sigma,  - (1:2)] <- calcCriteriaGR(contrast = contrast, groups = (1:n) %in% res.GR$GR, W = W, 
                                                         sigma = sigma[iter_sigma], breaks = breaks, rm.warning = FALSE, 
                                                         criterion.transition = criterion.transition, criterion.sdfront = criterion.sdfront, 
                                                         criterion.entropy = criterion.entropy, criterion.Kalinsky = criterion.Kalinsky, criterion.Laboure = criterion.Laboure)[,-1]
    
    ## save optimum
    for(iter_criteria in 1:n.criteria){
      if(!is.na(df.criterion[iter_sigma, iter_criteria + 2]) && df.criterion[iter_sigma, iter_criteria + 2] > best[iter_criteria]){
        list.GR[[iter_criteria]] <- res.GR$GR
        best[iter_criteria] <- df.criterion[iter_sigma, iter_criteria + 2]
      }            
    }
    
    ## update
    if(length(res.GR$GR) > 0){seed.sigma <- res.GR$GR}    
    if(length(res.GR$GR) == n){break}
  }
  if(verbose == TRUE){cat("\n")}
  
  return(list(df.criterion = df.criterion, 
              list.GR = list.GR, 
              best = best, 
              n.max = length(contrast))
  )
}

GRalgo <- function(contrast, W, seed, sigma_max, range, breaks, step, operator, iter_max, keep.lower, keep.upper, 
                   history.sigma, history.step, history.front){  
  
  #### initialisation 
  if(history.sigma){
    hist_sigma <- matrix(NA, nrow = iter_max, ncol = 2)
  }
  if(history.step){
    hist_stepGR <- rep(0, length(seed))
  }
  if(history.front){
    hist_frontGR <- rep(0, length(seed))
  }
  
  iter <- 0
  test.break <- FALSE  
  
  GR <- seed
  GRsave <- NULL
  breaks.GR <- breaks[graphics::hist(contrast[GR], breaks = breaks, plot = FALSE)$density > 0]
  n.breaks.GR <- length(breaks.GR)
  
  ##### definition of the GR algorithm
  
  fct_while <- "function(){"
  fct_while <- paste(fct_while, 
                     "while(identical(GR, GRsave) == FALSE && iter < iter_max){", 
                     "iter <- iter + 1", 
                     "GRsave <- GR \n", 
                     sep = "\n")
  
  ## dilatation
  fct_while <- paste(fct_while, "## dilatation", 
                     "Dnew <- unique(c(GR, which(spam::colSums(W[GR,, drop = FALSE]) > 0)))", 
                     "Dnew <- Dnew[ (contrast[Dnew] >= range[1]) * (contrast[Dnew] <= range[2]) == 1 ] # supplement perso pour ne pas ajouter des px en dehors du seuil possible", 
                     "if(length(Dnew) <= 1){GR <- Dnew ; test.break = TRUE ; break} # check sigma validity", 
                     "# sigma", 
                     paste("sigma <- ", operator, "(contrast[Dnew])", sep = ""), 
                     if(history.sigma == TRUE){"hist_sigma[iter,1] <- sigma \n"}, 
                     sep = "\n")
  
  ## homogeneous         
  fct_while <- paste(fct_while, "## homogeneous", 
                     "if(sigma <= sigma_max){", 
                     "GR <- Dnew", 
                     "}", 
                     sep = "\n")
  
  ## reduce dilatation (breaks correspond to Inf)
  fct_while <- paste(fct_while, 
                     "else{ ## reduce dilatation", 
                     sep = "")
  
  fct_while <- paste(fct_while, 
                     "Dnew_red <- Dnew[ (contrast[Dnew] >= breaks.GR[1]) * (contrast[Dnew] <= (breaks.GR[length(breaks.GR)] + step)) == 1 ]", 
                     "# sigma", 
                     "if(length(Dnew_red) <= 1){ test.break = TRUE ; break}", 
                     paste("sigma <- ", operator, "(contrast[Dnew_red]) \n", sep = ""), 
                     sep = "\n")
  
  ## homogenous : extensions
  fct_while <- paste(fct_while, "## homogenous : extensions", 
                     "if(sigma <= sigma_max){", 
                     "Neighborhood <- setdiff(which(spam::colSums(W[Dnew_red,, drop = FALSE]) > 0), Dnew_red)", 
                     if(keep.lower == TRUE || keep.upper == FALSE){"test.leftExtension <- (contrast[Neighborhood] >= (breaks.GR[1]-step)) * (contrast[Neighborhood] < (breaks.GR[n.breaks.GR] + step))"}, 
                     if(keep.lower == TRUE || keep.upper == FALSE){"leftExtension <- Neighborhood[which(test.leftExtension == 1)]"}, 
                     if(keep.upper == TRUE || keep.lower == FALSE){"test.rightExtension  <- (contrast[Neighborhood] >= (breaks.GR[1])) * (contrast[Neighborhood] < (breaks.GR[n.breaks.GR] + 2 * step))"}, 
                     if(keep.upper == TRUE || keep.lower == FALSE){"rightExtension <- Neighborhood[which(test.rightExtension == 1)]"}, 
                     "# sigma", 
                     sep = "\n")
  
  ## full extension
  if(keep.lower == keep.upper){
    fct_while <- paste(fct_while, "# full extension", 
                       paste("sigma <- ", operator, "(contrast[c(Dnew_red, leftExtension, rightExtension)])", sep = ""),                     
                       "if(sigma <= sigma_max){", 
                       "GR <- c(Dnew_red, unique(c(leftExtension, rightExtension)))", 
                       "} else {", 
                       sep = "\n")
  }
  
  ## left extension
  if(keep.lower == TRUE || keep.upper == FALSE){
    fct_while <- paste(fct_while, "# left extension", 
                       "if(length(leftExtension) > 0){", #  && breaks.GR[1]-step >= range[1]                     
                       paste("sigmaL <- ", operator, "(contrast[c(Dnew_red, leftExtension)])", sep = ""), 
                       "}else{sigmaL <- sigma_max + 1}", 
                       sep = "\n")
  }
  
  ## right extension
  if(keep.upper == TRUE || keep.lower == FALSE){
    fct_while <- paste(fct_while, "# right extension", 
                       "if(length(rightExtension) > 0){", #  breaks.GR[n.breaks.GR] + 2 * step <= range[2]            
                       paste("sigmaR <- ", operator, "(contrast[c(Dnew_red, rightExtension)])", sep = ""), 
                       "} else {sigmaR <- sigma_max + 1}", 
                       sep = "\n")
  }           
  
  if(keep.lower == TRUE || keep.upper == FALSE){
    fct_while <- paste(fct_while, 
                       paste("if(sigmaL <= sigma_max", if(keep.upper == keep.lower){" && sigmaL <= sigmaR"}, "){", sep = ""), 
                       "GR <- c(Dnew_red, leftExtension)", 
                       if(keep.upper == keep.lower){"} else "}else{"} else {"}, 
                       sep = "\n")
  }
  
  if(keep.upper == TRUE || keep.lower == FALSE){
    fct_while <- paste(fct_while, 
                       "if(sigmaR <= sigma_max){", 
                       "GR <- c(Dnew_red, rightExtension)", 
                       "} else {", 
                       sep = "\n")
  }
  
  fct_while <- paste(fct_while, 
                     "GR <- Dnew_red", 
                     "}", 
                     if(keep.lower == keep.upper){"}"}, 
                     "}", 
                     sep = "\n")
  
  ## contraction
  if(keep.lower == FALSE || keep.upper == FALSE){
    fct_while <- paste(fct_while, 
                       "else {", 
                       sep = "")
    
    fct_while <- paste(fct_while, "## contraction", 
                       "ordreC <- 1",      
                       "max <- 0 \n", 
                       "while(max == 0){", 
                       "if(ordreC > n.breaks.GR){GR <- GRsave ;  test.break = TRUE ; break} \n", 
                       "for(iter_C in 0:ordreC){", 
                       paste("breaks.min <- breaks.GR[1]", if(keep.lower == FALSE){"+ iter_C * step"}, sep = ""), 
                       paste("breaks.max <- breaks.GR[n.breaks.GR] + step", if(keep.upper == FALSE){" - (ordreC - iter_C) * step"}, sep = ""), # le + step est un ajout
                       "C <- Dnew_red[ (contrast[Dnew_red] >= breaks.min) * (contrast[Dnew_red] <= breaks.max) == 1 ]", 
                       "if(length(C) <= 1){GR <- GRsave ;  test.break = TRUE ; break }", 
                       paste("sigma.test <- ", operator, "(contrast[C])", sep = ""),            
                       "if(sigma.test <= sigma_max && length(C) > max){", 
                       "max <- length(C)", 
                       "GR <- C", 
                       "}", 
                       "}", 
                       "if(max == 0){ordreC <- ordreC + 1}", 
                       "} \n", 
                       "} # end non homogeneous dilatation", 
                       sep = "\n")
  }
  
  ## end loop
  fct_while <- paste(fct_while,                    
                     "} # end non homogeneous", 
                     "GR <- sort(GR)", 
                     "breaks.GR <- breaks[graphics::hist(contrast[GR], breaks = breaks, plot = FALSE)$density > 0]", 
                     "n.breaks.GR <- length(breaks.GR) \n \n", 
                     sep = "\n")
  
  ## save
  if(history.sigma == TRUE){
    fct_while <- paste(fct_while, 
                       "hist_sigma[iter, 2] <- stats::sd(contrast[GR])", 
                       sep = "\n")
  }
  
  if(history.step == TRUE || history.front == TRUE){
    fct_while <- paste(fct_while, 
                       "GR_in_GRsave <- which(GR %in% GRsave)", 
                       "GRsave_in_GR <- which(GRsave %in% GR)", 
                       sep = "\n")
    
    if(history.step == TRUE){
      fct_while <- paste(fct_while, 
                         "hist_stepGR.sauve <- hist_stepGR", 
                         "hist_stepGR <- rep(iter, length(GR))", 
                         "hist_stepGR[GR_in_GRsave] <- hist_stepGR.sauve[GRsave_in_GR]", 
                         sep = "\n")
    }
    
    if(history.front == TRUE){
      fct_while <- paste(fct_while, 
                         "hist_frontGR.sauve <- hist_frontGR", 
                         "hist_frontGR <- rep(-1, length(GR)) \n", 
                         "# old GR points", 
                         "hist_frontGR[GR_in_GRsave] <- hist_frontGR.sauve[GRsave_in_GR] \n", 
                         "# new GR points", 
                         "newGR <- setdiff(GR, GRsave)", 
                         "if(length(newGR) > 0){", 
                         "front <- as.factor(calcGroupsW_cpp(W_i = W@i, W_p = W@p, newGR - 1, 10000, false)$group_subset) \n", 
                         "# order front by ", 
                         "levels(front) <- levels(front)[rank(-table(front) + seq(0, 0.1, length.out = length(levels(front))))]", 
                         "front <- as.numeric(as.character(front))", 
                         "if(iter == 1){", 
                         "hist_frontGR[GR %in% GRsave == 0] <- front", 
                         "}else{", 
                         "origin <- sapply(newGR, function(x){  # privilegie l appartenance aux grands fronts i.e. max(front)", 
                         "index_Wn0 <- which(W[x, GRsave] > 0)", 
                         "valW <- W[x, GRsave[index_Wn0]]", 
                         "if(length(valW) > 0){max(hist_frontGR.sauve[index_Wn0[which(valW == min(valW))]])}else{\"\"}", 
                         "}) \n",            
                         "newGR2 <- newGR[which(nchar(origin) == 0)] \n",            
                         "if(length(newGR2) > 0){ # cas ou il y a eu une extension il faut recuperer le voisinage precedent", 
                         if(keep.upper == keep.lower){"newGR_noExp <- which(newGR %in% c(leftExtension, rightExtension) == FALSE)"}, 
                         if(keep.upper == TRUE && keep.lower == FALSE){"newGR_noExp <- which(newGR %in% rightExtension == FALSE)"}, 
                         if(keep.lower == TRUE && keep.upper == FALSE){"newGR_noExp <- which(newGR %in% leftExtension == FALSE)"}, 
                         "restrictedGR <- newGR[newGR_noExp] \n", 
                         "origin2 <- sapply(newGR2, function(x){", 
                         "index_Wn0 <- which(W[x,restrictedGR] > 0)", 
                         "W_red <- W[x, restrictedGR[index_Wn0]]", 
                         "max(origin[newGR_noExp[index_Wn0[which(W_red == min(W_red))]]])", 
                         "})", 
                         "origin[newGR %in% newGR2] <- origin2", 
                         "} \n",            
                         "hist_frontGR[GR %in% newGR] <- paste(origin, front, sep = \".\")", 
                         "}", 
                         "}", 
                         "hist_frontGR[GR_in_GRsave] <- paste(hist_frontGR[GR_in_GRsave], \".\", sep = \"\")", 
                         sep = "\n")
    } # end if history.front
  } # end if history.step || history.front
  
  ## export
  fct_while <- paste(fct_while, "## contraction", 
                     "} # end while loop \n", 
                     " ## export ", 
                     "return(list(GR = GR, ", 
                     "test.break = test.break, ", 
                     "iter = iter, ", 
                     "test.id = identical(GR, GRsave), ", 
                     if(history.sigma == TRUE){"sigma = if(iter > 0){hist_sigma[1:iter,]}else{NULL}, "}else{"sigma = NULL, "}, 
                     if(history.step == TRUE){"history.step = if(iter > 0){hist_stepGR}else{NULL}, "}else{"history.step = NULL, "}, 
                     if(history.front == TRUE){"history.front = if(iter > 0){hist_frontGR}else{NULL}"}else{"history.step = NULL"}, 
                     ")
                     ) \n",                     
                     "} # end function ", 
                     sep = "\n")
  
  while_loop <- eval(parse(text = paste(fct_while)))

  #### computation and export
  
  return(while_loop())
}

#### 2- plot function

plotSigmaGR <- function(calcSigmaGR, mar = c(3, 3, 2, 2), mgp = c(2, 0.75, 0), main = "", col = c("black", grDevices::rainbow(5)), 
                        criterion = c("n.obs", "transition", "sdfront", "entropy", "Kalinsky", "Laboure"), 
                        name_criteria = c("region size", "boundary transition", "boundary heterogeneity", 
                                          "region entropy", "region Kalinsky", "region Laboure"), 
                        filename = "calcSigmaGR", ...){
  
  #### pre initialisation 
  ## get graphical parameters 
  dots.arguments <- list(...)
  names_dots.arguments <- names(dots.arguments)
  
  ## get data  
  data <- calcSigmaGR$df.criterion
  
  #### tests
  valid.values <- c("n.obs","transition","sdfront","entropy","Kalinsky","Laboure") # to be checked !!!!
  
  if (optionsMRIaggr("checkArguments")) {
  
  validNumeric(value = mar, validLength = 4, method = "plotSigmaGR")
  validNumeric(value = mgp, validLength = 3, method = "plotSigmaGR")
  validCharacter(value = main, validLength = 1, method = "plotSigmaGR")
  validCharacter(value = col, validLength = 6, method = "plotSigmaGR")
  validCharacter(value = criterion, validLength = NULL, 
                 validValues = c("n.obs", "transition", "sdfront", "entropy", "Kalinsky", "Laboure"), 
                 method = "plotSigmaGR")
  validCharacter(value = name_criteria, validLength = 6, method = "plotSigmaGR")
  validCharacter(value = filename, validLength = 1, method = "plotSigmaGR")
  validCharacter(names_dots.arguments, name = "...", validLength = NULL, 
                 validValues = c("window", "width", "height", "path", "unit", "res"), 
                 refuse.NULL = FALSE, method = "plotSigmaGR")
  
  }
    
  #### initialisation
  ## set specific display
  optionsMRIaggr.eg <- optionsMRIaggr()[c("window", "width", "height", "path", "unit", "res")]        
  if(length(names_dots.arguments) > 0){
    optionsMRIaggr.eg[names_dots.arguments] <- dots.arguments[names_dots.arguments]
  }
  window <- optionsMRIaggr.eg$window
  width <- optionsMRIaggr.eg$width
  height <- optionsMRIaggr.eg$height
  path <- optionsMRIaggr.eg$path
  unit <- optionsMRIaggr.eg$unit
  res <- optionsMRIaggr.eg$res
  
  res.initW <- initWindow(window, filename, path, width, height, unit, res, 
                          n.plot = 1, mfrow = c(1, 1), xlim = NULL, ylim = NULL, method = "calcSigmaGR")
  scale <- res.initW$scale

  ## restriction to the requested criteria
  index_criteria <- which(names(data)[-1] %in% intersect(criterion,valid.values))
  subset_criteria <- names(data)[-1][index_criteria]
  n.criterion <- length(subset_criteria)
  name_criteria <- name_criteria[index_criteria]
  col <- col[index_criteria]
 
  df.criterion <- data[,subset_criteria, drop = FALSE]

  #### display
  if (all(is.na(data)[,-1])) {
    warning("plotSigmaGR : all criteria are NA, plot is not possible \n")
  }else{
    
    # normalisation
    if ("n.obs" %in% subset_criteria) {df.criterion$n.obs <- df.criterion$n.obs / calcSigmaGR$n.max}
    
    if (any(valid.values[-1] %in% subset_criteria)) {
      df.criterion[,setdiff(subset_criteria, "n.obs")] <- sweep(df.criterion[,setdiff(subset_criteria, "n.obs"), drop = FALSE], 2, FUN = "/", 
                                                                apply(df.criterion[,setdiff(subset_criteria, "n.obs"), drop = FALSE], 2, max, na.rm = TRUE))
    }
    
    
    initDisplayWindow(window = window, filename = filename, path = path, width = width, height = height, scale = scale, res = res, 
                      mfrow = c(1, 1), bg = graphics::par()$bg, pty = graphics::par()$pty, mar = mar, mgp = mgp)
    
    index_sigma <- which(rowSums(!is.na(df.criterion)) > 0)
    graphics::plot(x = data$sigma[index_sigma], y = rep(NA, length(index_sigma)), ylim = c(0, 1.25), 
                   xlab = expression(sigma[max]), ylab = "scaled value (% max)", main = main)
    
 
    for(iter_criterion in 1:n.criterion){      
      graphics::points(data$sigma, df.criterion[,iter_criterion], col = col[iter_criterion], 
                       type = "b", pch = 20)
      
      if(criterion[iter_criterion] != "n.obs"){
        optimum <- which.max(df.criterion[,iter_criterion])
        graphics::points(data$sigma[optimum], 1, col = col[iter_criterion], pch = 21, cex = 1.5)
        graphics::points(x = rep(data$sigma[optimum], 2), y = c(0,1), type = "l", col = col[iter_criterion], lty = 2)
      }
    }    
    
    graphics::legend("topleft", legend = name_criteria[1:min(n.criterion,3)], col = col[1:min(n.criterion,3)], pch = 20, bty = "n")
    if(n.criterion > 3){
      graphics::legend("topright", legend = name_criteria[4:min(n.criterion,6)], col = col[4:min(n.criterion,6)], pch = 20, bty = "n")
    }
    switch(window, 
           eps = grDevices::dev.off(), 
           svg = grDevices::dev.off(), 
           png = grDevices::dev.off()
    )
  }
}


#### 3- init functions ####
initGR <- function(contrast, W, seed, sigma_max, range, range.seed, breaks, iter_max, sd.robust,
                   keep.lower, keep.upper, history.sigma, history.step, history.front,
                   rescale,  verbose, method = "initGR"){
  
  #### pre intialization 
  if (is.logical(seed)) {seed <- which(seed)}  
  n <- length(contrast)
    
  #### tests 
  if (optionsMRIaggr("checkArguments")) {
    
  validNumeric(value = contrast, validLength = NULL, method = method)
  validDim_matrix(value1 = W, value2 = c(n,n), name1 = "W", method = method)
  validInteger(value = seed, validLength = NULL, min = 1, max = length(contrast), refuse.duplicates = TRUE, method = method)  
  validNumeric(value = sigma_max, validLength = 1, min = 0, method = method)
  validNumeric(value = range, validLength = 2, method = method)
  if( all( (contrast >= range[1]) * (contrast <= range[2]) == FALSE) ){
    stop(method, " : incorrect specification of \'range\' \n", 
         "requested contrast range : ", paste(range, collapse = " "), "\n", 
         "\'contrast\' range ", paste(range(contrast), collapse = " "), "\n")
  }
  validNumeric(value = range.seed, validLength = 2, method = method)
  test.seed <- (contrast[seed] >= range.seed[1]) * (contrast[seed] <= range.seed[2])
  if (all(test.seed == FALSE)) {
    stop(method, " : incorrect specification of \'range.seed\' \n", 
         "requested contrast range.seed : ", paste(range.seed, collapse = " "), "\n", 
         "contrast range of the \'seeds\' ", paste(range(contrast[seed]), collapse = " "), "\n")
  }
  validNumeric(value = breaks, validLength = NULL, method = method)
  validInteger(value = iter_max, validLength = 1, min = 0, method = method)
  validLogical(value = sd.robust, validLength = 1, method = method)
  validLogical(value = keep.lower, validLength = 1, method = method)
  validLogical(value = keep.upper, validLength = 1, method = method)
  validLogical(value = history.sigma, validLength = 1, method = method)
  validLogical(value = history.step, validLength = 1, method = method)
  validLogical(value = history.front, validLength = 1, method = method)
  validLogical(value = rescale, validLength = 1, method = method)
  validLogical(value = verbose, validLength = 1, method = method)

  }

  #### intialization
  if (rescale) {contrast <- as.numeric(scale(contrast))}
  
  if (length(breaks) == 1) {breaks <- seq(min(contrast), max(contrast), length.out = breaks)}
  step <- breaks[2] - breaks[1] + 10^{-12}
  
  seed <- seed[test.seed == 1]

  if (verbose == TRUE) {cat("number of valid seeds : ", sum(test.seed), " over ", length(seed), " seeds \n", sep = "")}
  
  #### export ####  
  return( list(contrast = contrast,
              breaks = breaks,
              step = step,
              seed = seed) )
}
