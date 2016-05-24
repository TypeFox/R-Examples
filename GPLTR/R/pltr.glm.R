pltr.glm <- function(data, Y.name, X.names, G.names, family = 'binomial', args.rpart = list(cp=0, minbucket=20, maxdepth=10), epsi = 1e-3, iterMax = 5, iterMin = 3, verbose = TRUE)
 { 
  time1 <- Sys.time()
  
  
  logistic.split <- list(eval = .logistic.eval, split = .logistic.split, init = .logistic.init)
  
  
  
  if(iterMin > iterMax) stop("iterMax not greater than iterMin !")
  
  ##	Method
  if(family == "gaussian") method = "anova"
  if(family == "binomial") method = logistic.split
  
  #	Initial fit 
  fit0 <- glm(as.formula(paste(Y.name, "~ -1+", paste(X.names, collapse=" + "))), data = data, family = family)
  
  hat_gamma <- fit0$coef
  
  #	Iterations
  diff_norm_gamma <- 10 # Valeur arbitraire
  
  ## Number of iterations
  nber_iter <- 1
  
  ## taking into account the matricial product
  if(length(X.names) > 1) product <- "%*%" 
  else product <- "*"
  
  if(verbose) cat("Iteration process...\n\n")
  
  while((diff_norm_gamma>0) & ((nber_iter<=iterMin) | (diff_norm_gamma>=epsi)) & (nber_iter<=iterMax))
  {
    #	1.a : fit tree
    fit_tree <- rpart(as.formula(paste(Y.name, " ~ ", paste(G.names, collapse=" + "), paste("+ offset(offsetX)"))), data = eval(parse(text = paste("data.frame(data, offsetX = as.matrix(data[,X.names])", product, "hat_gamma)"))), method = method, control = args.rpart)
    
    #	1.b : fit linear model
    indicators_tree <- sapply(tree2indicators(fit_tree), function(u) return(paste("as.integer(", u, ")")))
    
    nber_indicators <- length(indicators_tree)
    
    fit_lm <- glm(as.formula(paste(Y.name, "~ ", paste(indicators_tree[-nber_indicators], collapse = "+"), paste("+ offset(offsetX)"))), eval(parse(text = paste("data.frame(data, offsetX = as.matrix(data[,X.names])", product, "hat_gamma)"))), family = family)
    
    #	1.c : Update the estimate of gamma
    hat_beta <- fit_lm$coef
    
    offsetZ <- with(data, eval(parse(text = paste(hat_beta, c(1, indicators_tree[-nber_indicators]), collapse = " + ", sep="*")))) ##
    
    fit_lm_update <- glm(as.formula(paste(Y.name, "~ -1 + ", paste(X.names, collapse = " + "), paste("+ offset(offsetZ)"))), data = data.frame(data, offsetZ = offsetZ), family = family)
    
    hat_gamma_update <- fit_lm_update$coef
    
    #	2 : Stop conditions
    diff_norm_gamma <- .norm2(hat_gamma_update - hat_gamma)
    
    hat_gamma <- hat_gamma_update
    
    ##	Number of iterations
    if(verbose) cat("Iteration ", nber_iter, "in PLTR; Diff_norm_gamma = ", diff_norm_gamma, "\n")
    nber_iter <- nber_iter + 1
  }
  if(verbose)
    { 
     cat("End of iteration process\n")
     cat("Number of iterations: ", nber_iter - 1, "\n\n")
    }
  time2 <- Sys.time()
  Timediff <- difftime(time2, time1)
  return(list(fit = fit_lm_update, tree = fit_tree, nber_iter = nber_iter - 1, Timediff = Timediff))
}
