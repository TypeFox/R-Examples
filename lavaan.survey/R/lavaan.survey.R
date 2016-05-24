# Complex sampling analysis of SEM models
# Daniel Oberski, 2015-11-03

lavaan.survey <- 
  function(lavaan.fit, survey.design, 
           estimator=c("MLM", "MLMV", "MLMVS", "WLS", "DWLS", "ML"),
           estimator.gamma=c("default","Yuan-Bentler")) {
  
  # Not all estimators in lavaan make sense to use here, therefore matching args
  estimator <- match.arg(estimator) 
  if(estimator=="ML") warning("Estimator 'ML' will not correct standard errors and chi-square statistic.")
  estimator.gamma <- match.arg(estimator.gamma) # Smoothing Gamma or not
  
  # Names of the observed variables (same for each group)
  ov.names <- lavaanNames(lavaan.fit, type="ov", group=1)
  
  # The MP-inverse duplication matrix is handy for removing redundancy
  Dplus <- lavaan::lav_matrix_duplication_ginv( length(ov.names) )
  # Create a formula that includes all observed variables for svymean
  ov.formula <- as.formula(paste("~", paste(ov.names, collapse="+")))
  
  # <no. group>-sized lists that will contain the asy. covariance matrix,
  #  and sample covariance matrix and mean vector
  ngroups <- lavInspect(lavaan.fit, "ngroups")
  
  Gamma <- vector("list", ngroups)
  sample.cov <- vector("list", ngroups)
  sample.mean <- vector("list", ngroups)
  sample.nobs <- vector("list", ngroups)
  
  for(g in seq(ngroups)) {
    if(ngroups > 1) {
      # Use survey::subset to create data groups
      survey.design.g <- 
        subset(survey.design, eval(parse(text=sprintf("%s == '%s'", 
                                                      lavInspect(lavaan.fit, "group"), 
                                                      lavInspect(lavaan.fit, "group.label")[[g]]))))
    } 
    else { # In case of no groups, just use the original survey design object.
      survey.design.g <- survey.design  
    }
    
    # Function that takes survey design and returns the Gamma & observed moments
    get.stats.design <- function(survey.design.g, sample.nobs.g) {
      sample.cov.g <- as.matrix(svyvar(ov.formula, design=survey.design.g, na.rm=TRUE))  
      # survey package returns the variance matrix of the (co)variances as attr:
      Gamma.cov.g <- attr(sample.cov.g, "var")
      # Remove (co)variances wrt redundant elements of covs; not used by lavaan. 
      Gamma.cov.g <- Dplus %*% Gamma.cov.g %*% t(Dplus)
      
      # Same for mean vector
      sample.mean.g <- svymean(ov.formula, design=survey.design.g, na.rm=TRUE)  
      Gamma.mean.g <- attr(sample.mean.g, "var")
      
      # Join asy. variance matrices for means and covariances
      # TODO add offdiag
      Gamma.g <- lavaan::lav_matrix_bdiag(Gamma.mean.g, Gamma.cov.g)
      
      Gamma.g <- Gamma.g * sample.nobs.g # lavaan wants nobs * Gamma.
      
      # Since the above nonparametric estimate of Gamma can be unstable, Yuan
      # and Bentler suggested a model-smoothed estimate of it, optional here:
      if(estimator.gamma == "Yuan-Bentler") {
        r <- get.residuals(lavaan.fit) # Iff these asy = 0, all will be well...
        Gamma.g <- Gamma.g + (sample.nobs.g/(sample.nobs.g - 1)) * (r %*% t(r))
      }
      # Get rid of attribute, preventing errors caused by lazy evaluation
      # (This has to be at the end or lazy evaluation mayhem will ensue)
      attr(sample.cov.g, "var") <- NULL
      tmp  <- as.vector(sample.mean.g)
      names(tmp) <- names(sample.mean.g)
      sample.mean.g <- tmp	
      
      list(Gamma.g=Gamma.g, sample.cov.g=sample.cov.g, sample.mean.g=sample.mean.g)
    }
    # The data may be a list of multiply imputed datasets
    if(!any(class(survey.design.g) == "svyimputationList")) {
      # If no imputations, just use usual no. observations and asy variance
      sample.nobs.g <- lavInspect(lavaan.fit, "nobs")[[g]] 

      stats <- get.stats.design(survey.design.g, sample.nobs.g)
    } 
    else { # In case of multiply imputed data
      # Not only can nobs differ from lavaan.fit, but also per imputation
      sample.nobs.g <- get.sample.nobs(survey.design.g, lavInspect(lavaan.fit, "group"))
      
      # Retrieve point and variance estimates per imputation, 
      #    [TODO: this line will not be correct when nobs differs over imputations]
      stats.list <- lapply(survey.design.g[[1]], get.stats.design, sample.nobs=sample.nobs.g)
      m  <- length(stats.list) # no. imputation
      
      # Point estimates are average over imputations
      sample.cov.list <- lapply(stats.list, `[[`, 'sample.cov.g')
      sample.cov.g <- Reduce(`+`, sample.cov.list) / m
      cov.df <- Reduce(`rbind`, lapply(sample.cov.list, lavaan::lav_matrix_vech))
      sample.mean.list <- lapply(stats.list, `[[`, 'sample.mean.g')
      sample.mean.g <- Reduce(`+`, sample.mean.list) / m
      mean.df <- Reduce(`rbind`, sample.mean.list)
      
      # Variance estimates depend on within- and between-imputation variance:
      Gamma.within  <- Reduce(`+`, lapply(stats.list, `[[`, 'Gamma.g')) / m
      Gamma.between <- cov(cbind(mean.df, cov.df))
      Gamma.g <- Gamma.within + ((m + 1)/m) * Gamma.between
      
      # set stats with multiple imputation point and variance estimates
      stats <- list(Gamma.g=Gamma.g, sample.cov.g=sample.cov.g, sample.mean.g=sample.mean.g)
    }

    # Augment the list for this group
    Gamma[[g]] <- stats$Gamma.g
    sample.cov[[g]] <- stats$sample.cov.g
    sample.mean[[g]] <- stats$sample.mean.g
    sample.nobs[[g]] <- sample.nobs.g
  } # End of loop over groups

  new.call <- lavInspect(lavaan.fit, "call")
  new.call$data <- NULL                # Remove any data argument
  new.call$sample.cov <- sample.cov    # Set survey covariances
  new.call$sample.mean <- sample.mean  # Set survey means
  new.call$sample.nobs <- sample.nobs  
  new.call$estimator <- estimator  # Always use Satorra-Bentler or WLS estimator

  if(substr(estimator, 1, 2) == "ML") { # ML, robust options
    # Set asymptotic covariance matrix of sample means and covariances
    new.call$NACOV <- Gamma      
  }
  if(estimator %in% c("WLS", "DWLS")) {
    # Weighted Least Squares, adjust the weight matrix: MP inverse of Gamma
    # Note that Gamma may be singular.
    new.call$WLS.V <- lapply(Gamma, ginv)
  }
  new.fit <- eval(as.call(new.call)) # Run lavaan with the new arguments
  
  if(estimator %in% c("WLS", "DWLS")) return(new.fit) # We are done for WLS

  # For ML with robust se's, check that a possibly singular Gamma has not
  # created dependencies in the parameter estimates.
  # (Code below should really be implemented in lavaan...)
  evs.too.small <- sapply(Gamma, function(Gamma.g) {
    any(eigen(Gamma.g, only.values=TRUE)$values < (.Machine$double.eps*100))
  })
  if(any(evs.too.small)) {
    V.est <- lavaan::vcov(new.fit)
    if(any(Re(eigen(V.est, only.values=TRUE)$values) < (.Machine$double.eps*100))) {
      long.string  <- sprintf("Some of the standard errors may not be trustworthy.
        Some of the observed covariances or means are
        collinear, and this has generated collinearity in your
        parameter estimates.  This may be a sample size issue,
        missing data problem, or due to having too few
        clusters relative to the number of parameters. Problem
        encountered in group(s) %s",
        paste(which(evs.too.small), collapse=", "))
  
      warning(strwrap(long.string, width=9999, simplify=TRUE))#gotta love it
    }
  }

  new.fit
}

# Obtain residuals from a lavaan fit object, concatenating means w/ covariances
# (used in Yuan-Bentler correction)
get.residuals <- function(fit) {
    r  <- lavaan::residuals(fit)
    c(r$mean, lavaan::lav_matrix_vech(r$cov))
}

# Obtain sample size from multiply imputed svydesign object.
# In case sample size differs over imputations, takes median over imputations.
# TODO: Does not work with multiple group yet.
get.sample.nobs  <- function(svy.imp.design, group=NULL) {
  nobs.imp <- lapply(svy.imp.design[[1]], function(des) {nrow(des$variables)})
  return(median(unlist(nobs.imp)))
}

# Use the pFsum function from the survey package to obtain p value for the 
#   overall model fit using an F reference distribution where the 
#   denominator degrees of freedom is the design degrees of freedom.  
# An anonymous reviewer for J Stat Software suggested that 
#  "in surveys with small numbers of primary sampling units this sort of 
#   correction has often improved the 
#   behavior of tests in other contexts."
# The eigenvalues of the U.Gamma matrix will be the coefficients in the 
#   mixture of F's distribution (Skinner, Holt & Smith, pp. 86-87).
pval.pFsum <- function(lavaan.fit, survey.design, method = "saddlepoint") {
  # Check that Satorra-Bentler or Satterthwaite adjustment is present
  if(!lavInspect(lavaan.fit, "options")$test %in% 
              c("satorra.bentler", "mean.var.adjusted", "Satterthwaite")) {
    stop("Please refit the model with Satorra-Bentler (MLM) or Satterthwaite (MLMVS) adjustment.") 
  }
  
  UGamma <- lavTech(lavaan.fit, "ugamma")
  real.eigen.values <- Re(eigen(lavTech(lavaan.fit, "ugamma"), only.values = TRUE)$values)

  return(survey::pFsum(x=fitMeasures(lavaan.fit, "chisq"), df=rep(1, length(real.eigen.values)), 
                a=real.eigen.values, ddf=survey::degf(survey.design), lower.tail=FALSE,
                method=method))
}
