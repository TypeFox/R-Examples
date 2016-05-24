##generic boot.wt
boot.wt <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
                    sort = TRUE, nsim = 100, ...){
  cand.set <- formatCands(cand.set)
  UseMethod("boot.wt", cand.set)
}

##default to indicate when object class not supported
boot.wt.default <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
                            sort = TRUE, nsim = 100, ...) {
    stop("\nFunction not yet defined for this object class\n")
  }



##aov
boot.wt.AICaov.lm <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
                              sort = TRUE, nsim = 100, ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##extract data from first model
  data <- eval(cand.set[[1]]$call$data, environment(formula(cand.set[[1]])))

  ##create vector to store top model
  top <- character(nsim)
  
  for(i in 1:nsim) {
    ##resample data
    new.data <- data[sample(x = rownames(data), replace = TRUE), ]

    ##store results
    results <- lapply(X = cand.set, FUN = function(j) update(j, data = new.data))

    ##compute AIC table
    out <- aictab(cand.set = results, modnames = modnames, second.ord = second.ord,
                  nobs = nobs, sort = TRUE)

    ##determine highest-ranked model
    top[i] <- as.character(out$Modnames[1])
  }

  ##compute selection frequencies for each model
  rel.freqs <- table(top)/nsim

  
  ##check whether all models appear in table
  if(length(cand.set) != length(rel.freqs)) {
    ##assign 0 freqs for models never appearing at first rank
    all.freqs <- rep(0, times = length(modnames))
    names(all.freqs) <- modnames
    ##iterate over observed models
    for (k in 1:length(rel.freqs)){
      all.freqs[names(rel.freqs)[k]] <- rel.freqs[k]
    }
    rel.freqs <- all.freqs
  }
    

  ##original model selection
  orig.aic <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                     nobs = nobs, sort = FALSE)
  
  ##sort table according to model names
  orig.aic <- orig.aic[order(orig.aic$Modnames), ]

  ##sort relative frequencies according to model names
  rel.freqs <- rel.freqs[order(names(rel.freqs))] 

  ##rename column for relative frequencies
  names(orig.aic)[7] <- "PiWt"
  
  ##add column for pi_i
  orig.aic$PiWt<- rel.freqs

  ##reorder if required
  if (sort) {
    orig.aic <- orig.aic[order(orig.aic[, 3]), ]
    #orig.aic$Cum.PiWt <- cumsum(orig.aic[, 7])
  }
  
  class(orig.aic) <- c("boot.wt", "data.frame")
  return(orig.aic)
}



##betareg
boot.wt.AICbetareg <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
                               sort = TRUE, nsim = 100, ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##extract data from first model
  data <- eval(cand.set[[1]]$call$data, environment(formula(cand.set[[1]])))

  ##create vector to store top model
  top <- character(nsim)
  
  for(i in 1:nsim) {
    ##resample data
    new.data <- data[sample(x = rownames(data), replace = TRUE), ]

    ##store results
    results <- lapply(X = cand.set, FUN = function(j) update(j, data = new.data))

    ##compute AIC table
    out <- aictab(cand.set = results, modnames = modnames, second.ord = second.ord,
                  nobs = nobs, sort = TRUE)

    ##determine highest-ranked model
    top[i] <- as.character(out$Modnames[1])
  }

  ##compute selection frequencies for each model
  rel.freqs <- table(top)/nsim

  
  ##check whether all models appear in table
  if(length(cand.set) != length(rel.freqs)) {
    ##assign 0 freqs for models never appearing at first rank
    all.freqs <- rep(0, times = length(modnames))
    names(all.freqs) <- modnames
    ##iterate over observed models
    for (k in 1:length(rel.freqs)){
      all.freqs[names(rel.freqs)[k]] <- rel.freqs[k]
    }
    rel.freqs <- all.freqs
  }
    

  ##original model selection
  orig.aic <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                     nobs = nobs, sort = FALSE)

  ##sort table according to model names
  orig.aic <- orig.aic[order(orig.aic$Modnames), ]

  ##sort relative frequencies according to model names
  rel.freqs <- rel.freqs[order(names(rel.freqs))] 

  ##rename column for relative frequencies
  names(orig.aic)[7] <- "PiWt"
  
  ##add column for pi_i
  orig.aic$PiWt<- rel.freqs

  ##reorder if required
  if (sort) {
    orig.aic <- orig.aic[order(orig.aic[, 3]), ]
    #orig.aic$Cum.PiWt <- cumsum(orig.aic[, 7])
  }
  
  class(orig.aic) <- c("boot.wt", "data.frame")
  return(orig.aic)
}



##clm
boot.wt.AICsclm.clm <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
                               sort = TRUE, nsim = 100, ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##extract data from first model
  data <- eval(cand.set[[1]]$call$data, environment(formula(cand.set[[1]])))

  ##create vector to store top model
  top <- character(nsim)
  
  for(i in 1:nsim) {
    ##resample data
    new.data <- data[sample(x = rownames(data), replace = TRUE), ]

    ##store results
    results <- lapply(X = cand.set, FUN = function(j) update(j, data = new.data))

    ##compute AIC table
    out <- aictab(cand.set = results, modnames = modnames, second.ord = second.ord,
                  nobs = nobs, sort = TRUE)

    ##determine highest-ranked model
    top[i] <- as.character(out$Modnames[1])
  }

  ##compute selection frequencies for each model
  rel.freqs <- table(top)/nsim

  
  ##check whether all models appear in table
  if(length(cand.set) != length(rel.freqs)) {
    ##assign 0 freqs for models never appearing at first rank
    all.freqs <- rep(0, times = length(modnames))
    names(all.freqs) <- modnames
    ##iterate over observed models
    for (k in 1:length(rel.freqs)){
      all.freqs[names(rel.freqs)[k]] <- rel.freqs[k]
    }
    rel.freqs <- all.freqs
  }
    

  ##original model selection
  orig.aic <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                     nobs = nobs, sort = FALSE)

  ##sort table according to model names
  orig.aic <- orig.aic[order(orig.aic$Modnames), ]

  ##sort relative frequencies according to model names
  rel.freqs <- rel.freqs[order(names(rel.freqs))] 

  ##rename column for relative frequencies
  names(orig.aic)[7] <- "PiWt"
  
  ##add column for pi_i
  orig.aic$PiWt<- rel.freqs

  ##reorder if required
  if (sort) {
    orig.aic <- orig.aic[order(orig.aic[, 3]), ]
    #orig.aic$Cum.PiWt <- cumsum(orig.aic[, 7])
  }
  
  class(orig.aic) <- c("boot.wt", "data.frame")
  return(orig.aic)
}



##glm
boot.wt.AICglm.lm <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
                              sort = TRUE, nsim = 100, c.hat = 1, ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##extract data from first model
  data <- eval(cand.set[[1]]$call$data, environment(formula(cand.set[[1]])))

  ##create vector to store top model
  top <- character(nsim)
  
  for(i in 1:nsim) {
    ##resample data
    new.data <- data[sample(x = rownames(data), replace = TRUE), ]

    ##store results
    results <- lapply(X = cand.set, FUN = function(j) update(j, data = new.data))

    ##compute AIC table
    out <- aictab(cand.set = results, modnames = modnames, second.ord = second.ord,
                  nobs = nobs, sort = TRUE, c.hat = c.hat)

    ##determine highest-ranked model
    top[i] <- as.character(out$Modnames[1])
  }

  ##compute selection frequencies for each model
  rel.freqs <- table(top)/nsim

  
  ##check whether all models appear in table
  if(length(cand.set) != length(rel.freqs)) {
    ##assign 0 freqs for models never appearing at first rank
    all.freqs <- rep(0, times = length(modnames))
    names(all.freqs) <- modnames
    ##iterate over observed models
    for (k in 1:length(rel.freqs)){
      all.freqs[names(rel.freqs)[k]] <- rel.freqs[k]
    }
    rel.freqs <- all.freqs
  }
    

  ##original model selection
  orig.aic <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                     nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##sort table according to model names
  orig.aic <- orig.aic[order(orig.aic$Modnames), ]

  ##sort relative frequencies according to model names
  rel.freqs <- rel.freqs[order(names(rel.freqs))] 

  ##sort table according to model names
  orig.aic <- orig.aic[order(orig.aic$Modnames), ]

  ##sort relative frequencies according to model names
  rel.freqs <- rel.freqs[order(names(rel.freqs))] 

  ##rename column for relative frequencies
  names(orig.aic)[7] <- "PiWt"
  
  ##add column for pi_i
  orig.aic$PiWt<- rel.freqs

  ##reorder if required
  if (sort) {
    orig.aic <- orig.aic[order(orig.aic[, 3]), ]
    #orig.aic$Cum.PiWt <- cumsum(orig.aic[, 7])
  }
  
  class(orig.aic) <- c("boot.wt", "data.frame")
  return(orig.aic)
}



##hurdle models
boot.wt.AIChurdle <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
                                   sort = TRUE, nsim = 100, ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##extract data from first model
  data <- eval(cand.set[[1]]$call$data, environment(formula(cand.set[[1]])))

  ##create vector to store top model
  top <- character(nsim)
  
  for(i in 1:nsim) {
    ##resample data
    new.data <- data[sample(x = rownames(data), replace = TRUE), ]

    ##store results
    results <- lapply(X = cand.set, FUN = function(j) update(j, data = new.data))

    ##compute AIC table
    out <- aictab(cand.set = results, modnames = modnames, second.ord = second.ord,
                  nobs = nobs, sort = TRUE)

    ##determine highest-ranked model
    top[i] <- as.character(out$Modnames[1])
  }

  ##compute selection frequencies for each model
  rel.freqs <- table(top)/nsim

  
  ##check whether all models appear in table
  if(length(cand.set) != length(rel.freqs)) {
    ##assign 0 freqs for models never appearing at first rank
    all.freqs <- rep(0, times = length(modnames))
    names(all.freqs) <- modnames
    ##iterate over observed models
    for (k in 1:length(rel.freqs)){
      all.freqs[names(rel.freqs)[k]] <- rel.freqs[k]
    }
    rel.freqs <- all.freqs
  }
    

  ##original model selection
  orig.aic <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                     nobs = nobs, sort = FALSE)

  ##sort table according to model names
  orig.aic <- orig.aic[order(orig.aic$Modnames), ]

  ##sort relative frequencies according to model names
  rel.freqs <- rel.freqs[order(names(rel.freqs))] 

  ##rename column for relative frequencies
  names(orig.aic)[7] <- "PiWt"
  
  ##add column for pi_i
  orig.aic$PiWt<- rel.freqs

  ##reorder if required
  if (sort) {
    orig.aic <- orig.aic[order(orig.aic[, 3]), ]
    #orig.aic$Cum.PiWt <- cumsum(orig.aic[, 7])
  }
  
  class(orig.aic) <- c("boot.wt", "data.frame")
  return(orig.aic)
}



##lm
boot.wt.AIClm <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
                            sort = TRUE, nsim = 100, ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##extract data from first model
  data <- eval(cand.set[[1]]$call$data, environment(formula(cand.set[[1]])))

  ##create vector to store top model
  top <- character(nsim)
  
  for(i in 1:nsim) {
    ##resample data
    new.data <- data[sample(x = rownames(data), replace = TRUE), ]

    ##store results
    results <- lapply(X = cand.set, FUN = function(j) update(j, data = new.data))

    ##compute AIC table
    out <- aictab(cand.set = results, modnames = modnames, second.ord = second.ord,
                  nobs = nobs, sort = TRUE)

    ##determine highest-ranked model
    top[i] <- as.character(out$Modnames[1])
  }

  ##compute selection frequencies for each model
  rel.freqs <- table(top)/nsim

  
  ##check whether all models appear in table
  if(length(cand.set) != length(rel.freqs)) {
    ##assign 0 freqs for models never appearing at first rank
    all.freqs <- rep(0, times = length(modnames))
    names(all.freqs) <- modnames
    ##iterate over observed models
    for (k in 1:length(rel.freqs)){
      all.freqs[names(rel.freqs)[k]] <- rel.freqs[k]
    }
    rel.freqs <- all.freqs
  }
    

  ##original model selection
  orig.aic <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                     nobs = nobs, sort = FALSE)

  ##sort table according to model names
  orig.aic <- orig.aic[order(orig.aic$Modnames), ]

  ##sort relative frequencies according to model names
  rel.freqs <- rel.freqs[order(names(rel.freqs))] 

  ##rename column for relative frequencies
  names(orig.aic)[7] <- "PiWt"
  
  ##add column for pi_i
  orig.aic$PiWt<- rel.freqs

  ##reorder if required
  if (sort) {
    orig.aic <- orig.aic[order(orig.aic[, 3]), ]
    #orig.aic$Cum.PiWt <- cumsum(orig.aic[, 7])
  }
  
  class(orig.aic) <- c("boot.wt", "data.frame")
  return(orig.aic)
}



##multinom
boot.wt.AICmultinom.nnet <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
                                        sort = TRUE, nsim = 100, c.hat = 1, ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##extract data from first model
  data <- eval(cand.set[[1]]$call$data, environment(formula(cand.set[[1]])))

  ##create vector to store top model
  top <- character(nsim)
  
  for(i in 1:nsim) {
    ##resample data
    new.data <- data[sample(x = rownames(data), replace = TRUE), ]

    ##store results
    results <- lapply(X = cand.set, FUN = function(j) update(j, data = new.data))

    ##compute AIC table
    out <- aictab(cand.set = results, modnames = modnames, second.ord = second.ord,
                  nobs = nobs, sort = TRUE, c.hat = c.hat)

    ##determine highest-ranked model
    top[i] <- as.character(out$Modnames[1])
  }

  ##compute selection frequencies for each model
  rel.freqs <- table(top)/nsim

  
  ##check whether all models appear in table
  if(length(cand.set) != length(rel.freqs)) {
    ##assign 0 freqs for models never appearing at first rank
    all.freqs <- rep(0, times = length(modnames))
    names(all.freqs) <- modnames
    ##iterate over observed models
    for (k in 1:length(rel.freqs)){
      all.freqs[names(rel.freqs)[k]] <- rel.freqs[k]
    }
    rel.freqs <- all.freqs
  }
    

  ##original model selection
  orig.aic <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                     nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##sort table according to model names
  orig.aic <- orig.aic[order(orig.aic$Modnames), ]

  ##sort relative frequencies according to model names
  rel.freqs <- rel.freqs[order(names(rel.freqs))] 

  ##rename column for relative frequencies
  names(orig.aic)[7] <- "PiWt"
  
  ##add column for pi_i
  orig.aic$PiWt<- rel.freqs

  ##reorder if required
  if (sort) {
    orig.aic <- orig.aic[order(orig.aic[, 3]), ]
    #orig.aic$Cum.PiWt <- cumsum(orig.aic[, 7])
  }
  
  class(orig.aic) <- c("boot.wt", "data.frame")
  return(orig.aic)
}



##polr
boot.wt.AICpolr <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
                               sort = TRUE, nsim = 100, ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##extract data from first model
  data <- eval(cand.set[[1]]$call$data, environment(formula(cand.set[[1]])))

  ##create vector to store top model
  top <- character(nsim)
  
  for(i in 1:nsim) {
    ##resample data
    new.data <- data[sample(x = rownames(data), replace = TRUE), ]

    ##store results
    results <- lapply(X = cand.set, FUN = function(j) update(j, data = new.data))

    ##compute AIC table
    out <- aictab(cand.set = results, modnames = modnames, second.ord = second.ord,
                  nobs = nobs, sort = TRUE)

    ##determine highest-ranked model
    top[i] <- as.character(out$Modnames[1])
  }

  ##compute selection frequencies for each model
  rel.freqs <- table(top)/nsim

  
  ##check whether all models appear in table
  if(length(cand.set) != length(rel.freqs)) {
    ##assign 0 freqs for models never appearing at first rank
    all.freqs <- rep(0, times = length(modnames))
    names(all.freqs) <- modnames
    ##iterate over observed models
    for (k in 1:length(rel.freqs)){
      all.freqs[names(rel.freqs)[k]] <- rel.freqs[k]
    }
    rel.freqs <- all.freqs
  }
    

  ##original model selection
  orig.aic <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                     nobs = nobs, sort = FALSE)

  ##sort table according to model names
  orig.aic <- orig.aic[order(orig.aic$Modnames), ]

  ##sort relative frequencies according to model names
  rel.freqs <- rel.freqs[order(names(rel.freqs))] 

  ##rename column for relative frequencies
  names(orig.aic)[7] <- "PiWt"
  
  ##add column for pi_i
  orig.aic$PiWt<- rel.freqs

  ##reorder if required
  if (sort) {
    orig.aic <- orig.aic[order(orig.aic[, 3]), ]
    #orig.aic$Cum.PiWt <- cumsum(orig.aic[, 7])
  }
  
  class(orig.aic) <- c("boot.wt", "data.frame")
  return(orig.aic)
}



##rlm
boot.wt.AICrlm.lm <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
                                 sort = TRUE, nsim = 100, ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##extract data from first model
  data <- eval(cand.set[[1]]$call$data, environment(formula(cand.set[[1]])))

  ##create vector to store top model
  top <- character(nsim)
  
  for(i in 1:nsim) {
    ##resample data
    new.data <- data[sample(x = rownames(data), replace = TRUE), ]

    ##store results
    results <- lapply(X = cand.set, FUN = function(j) update(j, data = new.data))

    ##compute AIC table
    out <- aictab(cand.set = results, modnames = modnames, second.ord = second.ord,
                  nobs = nobs, sort = TRUE)

    ##determine highest-ranked model
    top[i] <- as.character(out$Modnames[1])
  }

  ##compute selection frequencies for each model
  rel.freqs <- table(top)/nsim

  
  ##check whether all models appear in table
  if(length(cand.set) != length(rel.freqs)) {
    ##assign 0 freqs for models never appearing at first rank
    all.freqs <- rep(0, times = length(modnames))
    names(all.freqs) <- modnames
    ##iterate over observed models
    for (k in 1:length(rel.freqs)){
      all.freqs[names(rel.freqs)[k]] <- rel.freqs[k]
    }
    rel.freqs <- all.freqs
  }
    

  ##original model selection
  orig.aic <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                     nobs = nobs, sort = FALSE)

  ##sort table according to model names
  orig.aic <- orig.aic[order(orig.aic$Modnames), ]

  ##sort relative frequencies according to model names
  rel.freqs <- rel.freqs[order(names(rel.freqs))] 

  ##rename column for relative frequencies
  names(orig.aic)[7] <- "PiWt"
  
  ##add column for pi_i
  orig.aic$PiWt<- rel.freqs

  ##reorder if required
  if (sort) {
    orig.aic <- orig.aic[order(orig.aic[, 3]), ]
    #orig.aic$Cum.PiWt <- cumsum(orig.aic[, 7])
  }
  
  class(orig.aic) <- c("boot.wt", "data.frame")
  return(orig.aic)
}



##survreg
boot.wt.AICsurvreg <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
                              sort = TRUE, nsim = 100, ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##extract data from first model
  data <- eval(cand.set[[1]]$call$data, environment(formula(cand.set[[1]])))

  ##create vector to store top model
  top <- character(nsim)
  
  for(i in 1:nsim) {
    ##resample data
    new.data <- data[sample(x = rownames(data), replace = TRUE), ]

    ##store results
    results <- lapply(X = cand.set, FUN = function(j) update(j, data = new.data))

    ##compute AIC table
    out <- aictab(cand.set = results, modnames = modnames, second.ord = second.ord,
                  nobs = nobs, sort = TRUE)

    ##determine highest-ranked model
    top[i] <- as.character(out$Modnames[1])
  }

  ##compute selection frequencies for each model
  rel.freqs <- table(top)/nsim

  
  ##check whether all models appear in table
  if(length(cand.set) != length(rel.freqs)) {
    ##assign 0 freqs for models never appearing at first rank
    all.freqs <- rep(0, times = length(modnames))
    names(all.freqs) <- modnames
    ##iterate over observed models
    for (k in 1:length(rel.freqs)){
      all.freqs[names(rel.freqs)[k]] <- rel.freqs[k]
    }
    rel.freqs <- all.freqs
  }
    

  ##original model selection
  orig.aic <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                     nobs = nobs, sort = FALSE)

  ##sort table according to model names
  orig.aic <- orig.aic[order(orig.aic$Modnames), ]

  ##sort relative frequencies according to model names
  rel.freqs <- rel.freqs[order(names(rel.freqs))] 

  ##rename column for relative frequencies
  names(orig.aic)[7] <- "PiWt"
  
  ##add column for pi_i
  orig.aic$PiWt<- rel.freqs

  ##reorder if required
  if (sort) {
    orig.aic <- orig.aic[order(orig.aic[, 3]), ]
    #orig.aic$Cum.PiWt <- cumsum(orig.aic[, 7])
  }
  
  class(orig.aic) <- c("boot.wt", "data.frame")
  return(orig.aic)
}



##vglm
boot.wt.AICvglm <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
                               sort = TRUE, nsim = 100, c.hat = 1, ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##extract data from first model
  data <- eval(cand.set[[1]]@call$data, environment(formula(cand.set[[1]])))

  ##create vector to store top model
  top <- character(nsim)
  
  for(i in 1:nsim) {
    ##resample data
    new.data <- data[sample(x = rownames(data), replace = TRUE), ]

    ##store results
    results <- lapply(X = cand.set, FUN = function(j) update(j, data = new.data))

    ##compute AIC table
    out <- aictab(cand.set = results, modnames = modnames, second.ord = second.ord,
                  nobs = nobs, sort = TRUE, c.hat = c.hat)

    ##determine highest-ranked model
    top[i] <- as.character(out$Modnames[1])
  }

  ##compute selection frequencies for each model
  rel.freqs <- table(top)/nsim

  
  ##check whether all models appear in table
  if(length(cand.set) != length(rel.freqs)) {
    ##assign 0 freqs for models never appearing at first rank
    all.freqs <- rep(0, times = length(modnames))
    names(all.freqs) <- modnames
    ##iterate over observed models
    for (k in 1:length(rel.freqs)){
      all.freqs[names(rel.freqs)[k]] <- rel.freqs[k]
    }
    rel.freqs <- all.freqs
  }
    

  ##original model selection
  orig.aic <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                     nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##sort table according to model names
  orig.aic <- orig.aic[order(orig.aic$Modnames), ]

  ##sort relative frequencies according to model names
  rel.freqs <- rel.freqs[order(names(rel.freqs))] 

  ##rename column for relative frequencies
  names(orig.aic)[7] <- "PiWt"
  
  ##add column for pi_i
  orig.aic$PiWt<- rel.freqs

  ##reorder if required
  if (sort) {
    orig.aic <- orig.aic[order(orig.aic[, 3]), ]
    #orig.aic$Cum.PiWt <- cumsum(orig.aic[, 7])
  }
  
  class(orig.aic) <- c("boot.wt", "data.frame")
  return(orig.aic)
}



##zeroinfl
boot.wt.AICzeroinfl <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
                                   sort = TRUE, nsim = 100, ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##extract data from first model
  data <- eval(cand.set[[1]]$call$data, environment(formula(cand.set[[1]])))

  ##create vector to store top model
  top <- character(nsim)
  
  for(i in 1:nsim) {
    ##resample data
    new.data <- data[sample(x = rownames(data), replace = TRUE), ]

    ##store results
    results <- lapply(X = cand.set, FUN = function(j) update(j, data = new.data))

    ##compute AIC table
    out <- aictab(cand.set = results, modnames = modnames, second.ord = second.ord,
                  nobs = nobs, sort = TRUE)

    ##determine highest-ranked model
    top[i] <- as.character(out$Modnames[1])
  }

  ##compute selection frequencies for each model
  rel.freqs <- table(top)/nsim

  
  ##check whether all models appear in table
  if(length(cand.set) != length(rel.freqs)) {
    ##assign 0 freqs for models never appearing at first rank
    all.freqs <- rep(0, times = length(modnames))
    names(all.freqs) <- modnames
    ##iterate over observed models
    for (k in 1:length(rel.freqs)){
      all.freqs[names(rel.freqs)[k]] <- rel.freqs[k]
    }
    rel.freqs <- all.freqs
  }
    

  ##original model selection
  orig.aic <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                     nobs = nobs, sort = FALSE)

  ##sort table according to model names
  orig.aic <- orig.aic[order(orig.aic$Modnames), ]

  ##sort relative frequencies according to model names
  rel.freqs <- rel.freqs[order(names(rel.freqs))] 

  ##rename column for relative frequencies
  names(orig.aic)[7] <- "PiWt"
  
  ##add column for pi_i
  orig.aic$PiWt<- rel.freqs

  ##reorder if required
  if (sort) {
    orig.aic <- orig.aic[order(orig.aic[, 3]), ]
    #orig.aic$Cum.PiWt <- cumsum(orig.aic[, 7])
  }
  
  class(orig.aic) <- c("boot.wt", "data.frame")
  return(orig.aic)
}



print.boot.wt <-
  function(x, digits = 2, ...) {
    cat("\nModel selection based on ", colnames(x)[3], ":\n", sep = "")
    if (any(names(x) == "c_hat")) {cat("(c-hat estimate = ", x$c_hat[1], ")\n", sep = "")}
    cat("\n")

    ##check if Cum.Wt should be printed
    nice.tab <- cbind(x[, 2], x[, 3], x[, 4], x[, 6], x[, 7])
    colnames(nice.tab) <- c(colnames(x)[c(2, 3, 4, 6, 7)])
    rownames(nice.tab) <- x[, 1]
    print(round(nice.tab, digits = digits)) #select rounding off with digits argument
    cat("\n")
  }

