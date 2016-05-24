lmf <-
function(formula,
                age,
                year,
                data,
                na.action = na.exclude,
                method = c("BFGS"),
                control = list(maxit = 500, reltol = sqrt(.Machine$double.eps)),
                ...)
{
  # Start; Set startup time
  ret <- list(running.time = proc.time()[3])
  # End
  
  # Start; Create the element which later takes the optim running time
  ret$optim.time <- NA
  # End
  
  # Start; COllect the function call
  ret$call <- match.call()
  # End
  
  # Start; Name data "x"
  # Data
  x <- data
  # End
  
  # Start; Check for missing arguments and match arguments to formal arguments
  # in the call, then match variable names to formal variable names in
  # the data set.
  # Collect the function call not expanding dots (...)
  mf <- match.call(expand.dots = FALSE)
  #  Check for missing arguments
  needed <- c("formula", "age", "year", "data")
  if(!all(needed %in% names(mf)))
  {
    nam <- needed[!(needed %in% names(mf))]
    msg <- sprintf(ngettext(length(nam),
                            "Argument %s is missing with no default",
                            "Arguments %s are missing with no defaults"),
                   paste(sQuote(nam), collapse = ", "))
    stop(msg)
  }
  #  Match call with formal arguments
  m <- match(c("formula", "age", "year", "data"), names(mf), 0L)
  #  Order mf according to the matching
  mf <- mf[c(1L, m)]
  # Extract relevant terms from call
  alv <- all.vars(mf$formula)
  #Change names in data set
  names(x)[names(x) %in% alv[1]] <- "recruits"
  names(x)[names(x) %in% alv[2]] <- "survival"
  if(is.null(mf$age)){
    x$age <- rep(0, length(x[, 1])) 
  } else{
    if(any(names(x) %in% as.character(mf$age))){
      names(x)[names(x) %in% as.character(mf$age)] <- "age"              
    } else{
      stop("Invalid 'age' argument.")
    }
  }
  if(is.null(mf$year)){
    stop("You need to specify the name of the 'year' variable.")
  } else{
    if(any(names(x) %in% as.character(mf$year))){
      names(x)[names(x) %in% as.character(mf$year)] <- "year"
    } else{
      stop("Invalid 'year' argument.")
    }
  }
  # End
  
  # Start; Change the formula to have Wj, estimated below, as response
  formula <- update.formula(old = formula, new = Wj~.)
  # End
  
  # Start; Preliminary checks 1
  # Check the data set for "NA" in age or year and remove these
  if(any(is.na(x[, 2:3])))
  {
    x <- x[rowSums(is.na(x[, 2:3])) == 0, ]
    warning("Entities with 'NA' for 'year' and/or 'age' were excluded")
  }
  # Numbers of parameters in the linear model
  ret$npar <- length(attr(terms(formula), "term.labels")) + 1
  # Number of observations for each age and year combination
  x.year.age <- xtabs(~ x$year + x$age)
  # Remove age classes with less than two years of enough data
  xage <- colSums(x.year.age >= ret$npar)
  if (any(xage < 2))
  {
    xage <- names(xage[xage < 2])
    msg <- sprintf(ngettext(length(xage),
                            "The age class %s had too few observations and was excluded",
                            "The age classes %s had too few observations and were excluded"),
                   paste(sQuote(xage), collapse = ", "))
    warning(msg, domain = NA)
    x <- x[!x$age %in% names(xage[xage < ret$npar]), ]
  }
  # Check that there are at least two years with enough data on all
  # age classes to proceed
  if(sum(rowSums(x.year.age <= ret$npar) == 0) < 2)
    stop("There are too few years with enough data from all age classes")
  # Age classes in the data
  ret$uage <- sort(unique(x$age))
  # Number of age classes in the data
  ret$nage <-  length(ret$uage)
  # The final age class
  ret$maxage <- ret$uage[ret$nage]
  # End
  
  # Start; Estimate projection matrix, and eigen value and vectors.
  # Point estimates
  ret$l <- promat(procomp(x, uage = ret$uage), nage = ret$nage)
  ret <- c(ret, eigenl(ret$l))
  # End
  
  # Start; Preliminary checks 2
  # Check that there are enough observations for each combination of age
  # and year to estimate parameters. Exclude as needed
  if(any(x.year.age <= ret$npar))
  {
    xyear <- names(c(rowSums(x.year.age <= ret$npar) > 0)[
      rowSums(x.year.age <= ret$npar) > 0])
    msg <- sprintf(ngettext(length(xyear),
                            "The year %s had too few observations for at least one age class and was excluded",
                            "The years %s had too few observations for at least one age class and were excluded"),
                   paste(sQuote(xyear), collapse = ", "))
    warning(msg, domain = NA)
    x <- x[!(x$year %in% xyear), ]
  }
  # Years in the data
  ret$uyear <- unique(x$year)
  # Number of years with data
  ret$nyear <-  length(ret$uyear)
  # The number of observation in the data set
  ret$nobs <- dim(x)[1]
  # The number of observation in the data set for each age class
  ret$nobs.age <- colSums(x.year.age)
  # Series of numbers as long as the number of observations
  ret$indnr <- 1 : ret$nobs
  # End
  
  # Start; Estimate Wj for all individuals
  x$Wj <- rep(NA, ret$nobs)
  for (i in ret$indnr)
  {
    x$Wj[i] <- ret$v[1] * x$recruits[i] + ret$v[ifelse(
      x$age[i] == ret$maxage, which(x$age[i] == ret$uage),
      which(x$age[i] == ret$uage) + 1)] * x$survival[i]
  }
  # End
  
  # Start; For each year within age: estimate alphas, alpha covariance matrises,
  # residual standard errors, degrees of freedom and residuals
  # Construct objects with appropriate numbers of predefined positions
  ret$ajt <- lapply(1 : ret$nage, function(i) vector("list", ret$nyear))
  ret$Ajt.us <- lapply(1 : ret$nage, function(i) vector("list", ret$nyear))
  ret$sigma.djt <- lapply(1 : ret$nage, function(i) vector("numeric", ret$nyear))
  ret$dof <- lapply(1 : ret$nage, function(i) vector("numeric", ret$nyear))
  ret$res <- lapply(1 : ret$nage, function(i) vector("list", ret$nyear))
  ret$fit <- lapply(1 : ret$nage, function(i) vector("list", ret$nyear))
  ret$leverage <- lapply(1 : ret$nage, function(i) vector("list", ret$nyear))
  ret$cook <- lapply(1 : ret$nage, function(i) vector("list", ret$nyear))
  # Estimate ajt, Ajt.us, rse and dof within each age for all years seperatedly
  for (i in 1 : ret$nage)
  {
    for(j in 1 : ret$nyear)
    {
      lx <- lm.extract(formula = formula,
                       data = x[x$age == ret$uage[i] & x$year == ret$uyear[j], ],
                       na.action = na.action)
      # Collect ajt
      ret$ajt[[i]][[j]] <-  lx$ajt
      # Collect Ajt.us
      ret$Ajt.us[[i]][[j]] <-  lx$Ajt.us
      # Collect rse
      ret$sigma.djt[[i]][j] <-  lx$sigma.djt
      # Collect dof
      ret$dof[[i]][j] <-  lx$dof
      # Return residuals to data set
      ret$res[[i]][[j]] <- lx$res
      # Return fitted values to data set
      ret$fit[[i]][[j]] <- lx$fit
      # Return leverage to data set
      ret$leverage[[i]][[j]] <- lx$leverage
      # Return cooks distance to data set
      ret$cook[[i]][[j]] <- lx$cook
    }
  }
  # End
  
  # Start; Estimate the age specific components of the demographic variance
  # and the total demographic variance
  # Estimate the age specific demographic variance components
  ret$sigma2.dj <- mapply(function(a, b){sum(a^2 * b) / sum(b)},
                          ret$sigma.djt, ret$dof, SIMPLIFY = FALSE)
  # Calculate the degrees of freedom for each sigma2.dj
  ret$sigma2.dj.dof <- mapply(sum, ret$dof, SIMPLIFY = FALSE)
  # Calculate the sample size for each age class
  Nj <- mapply(function(a) sum(a + ret$npar), ret$dof, SIMPLIFY = FALSE)
  # Estimate the standard deviation of the age specific demographic
  # variance components
  ret$sigma2.dj.sd <- mapply(function(a, b, c, d){sqrt(a^2 / (a - 1) * (mean(
    unlist(b)^4) - mean(unlist(b)^2)^2) / (a - c * d)^2)},
                             Nj, ret$res, ret$npar, ret$nyear, SIMPLIFY = FALSE)
  # Estimate the total demographic variance by summing and weighing by the
  # stable age distribution
  ret$sigma2.d <- sum(ret$u * unlist(ret$sigma2.dj))
  # Calculate the degrees of freedom for sigma2.d
  ret$sigma2.d.dof <- sum(unlist(ret$sigma2.dj.dof))
  # Estimate the standard deviation of the total demographic variance
  ret$sigma2.d.sd <- sqrt(sum(ret$u^2 * unlist(ret$sigma2.dj.sd)^2))
  # End
  
  # Start; Scale the covariance matrices by multiplying with the age specific demographic
  # variance components
  # Scale the (unscaled) covariances by the correct
  # age specific demographic variance
  ret$Ajt <- mapply(function(a, b){lapply(a, function(x) x * b)},
                    ret$Ajt.us, ret$sigma2.dj, SIMPLIFY = FALSE)
  # End
  
  # Start; Combine age specific alphas for each trait within years by taking
  # the sum of the alphas weighted by the stable age distribution to obtain
  # single yearly estimates of alphas
  ret$at <- lapply(1 : ret$nyear, function(i)
  {colSums(t(sapply(ret$ajt, "[[", i)) * ret$u)})
  # End
  
  # Start; Combine the age specific covariance matrices for the betas within
  # years by weighting by the stable age distribution to obtain single yearly
  # estimates of these matrices
  ret$At <- lapply(1 : ret$nyear, function(i)
  {Reduce('+', mapply('*', sapply(ret$Ajt, '[', i), as.list(ret$u^2), SIMPLIFY = FALSE))})
  # End
  
  # Start; Estimation of fluctuating selection (i.e. estimate M, the
  # temporal covariance matrix)
  ret$optim.time <- proc.time()[3]
  ret <- c(ret, fs(At = ret$At, at = ret$at, npar = ret$npar,
                   nyear = ret$nyear, method = method,
                   control = control, ...))
  ret$optim.time <- c(proc.time()[3] - ret$optim.time)
  # End
  
  # Start; Estimate E(at|?t), the yearly alpha estimates corrected for sampling error
  if(all(!is.na(ret$M)))
  {
    ret$atC <- atCfn(aM = ret$aM, M = ret$M, At = ret$At, at = ret$at)
  }
  else
    ret$atC <- rep(NA, ret$npar)
  # End
  
  # Start; Estimate alphas and alphas covariance matrix under the assumption
  # of no fluctuating selection (M = 0)
  ret <- c(ret, nfs(At = ret$At, at = ret$at, npar = ret$npar, nyear = ret$nyear))
  # End
  
  # Start; Estimate the total environmenal variance in the population
  ret$sigma2.e <- var(sapply(ret$at, '[', 1))
  # End
  
  # Start; Calculate total running time
  ret$running.time <- proc.time()[3] - ret$running.time
  # End
  
  # Start; Attach the data set to the output
  ret$data.set <- x
  # End
  
  # Start; Return the results of the analyses to the user
  class(ret) <- "lmf"
  ret
  # End
}
