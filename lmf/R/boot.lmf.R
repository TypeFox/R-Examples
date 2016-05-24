boot.lmf <-
function(object,
                     nboot = 1000,
                     what = c("projection", "alpha", "H0", "all"),
                     asim = c("ordinary", "parametric"),
                     sig.dj = TRUE,
                     H0exp = list(alpha = NULL, M = NULL),
                     H0con = c("fs", "nfs", "ds", "nds"),
                     method = c("BFGS"),
                     control = list(maxit = 500, reltol = sqrt(.Machine$double.eps)),
                     ...)
{
  #Arguments:
  #object: a fitted model object
  #nboot: The numbers of bootstraps desired
  #what: What set of parameters to estimate. "projection" bootstraps
  #the projection matrix (l), lambda, stable age distribution (u) and
  #the reproductive values (v). "alpha" bootstraps demographic and
  #environmental variances, yearly covariance matrices (At) (if sig.dj = TRUE),
  #yearly alpha estimates (at), temporal covariance matrix (M), temporal
  #alpha estimates (a(M)), yearly alpha estimates corrected for sampling
  #error (atC), temporal covariance matrix assuming no fluctuation
  #selection (?t(M=0)) and temporal alpha estimates assuming no
  #fluctuating selection (a(M=0))
  #asim: The type of bootstrap method for the yearly alpha estimates
  #sig.dj: Include uncertainty in the estimation of the demographic
  #variance when bootstrapping alpha estimates? (TRUE/FALSE)
  #H0exp: A list with the first element containing the alpha values and
  #the second element containing the M matrix elements
  #under the null hypothesis
  #H0con: The conditions under which the null hypothesis should be
  #tested. Options are: "fs" = fluctuating selection,
  #"nfs" = no fluctuating selection, "ds" = directional
  #selection and "nds" = no directional selection. "nds" is not
  #implemented due to increased risk of Type I error
  #method: The type of algorithm to use for the likelihood maximation.
  #Options are: "Nelder-Mead", "BFGS", "CG", "L-BFGS-B"
  #and "SANN". See ?optim for details
  #control = A list of control parameters for the maximization of the likelihood.
  #'maxit' controls the maximum number of iterations to use before convergence,
  #'reltol' controls the relative threshold for improvement in the likelihood which
  #desides whether to continue maximation or end. see ?optim for details
  #...: Alternative options passed to optim
  #Set the starting running time
  ret <- list(running.time = proc.time()[3])
  #Keep the call
  ret$call <- match.call()
  #Rename model output
  x <- object
  #Check model class
  if(class(x) != "lmf")
    stop("The model is not of class 'lmf'")
  #Choose "all" if all or no options for what to bootstrap are choosen
  ifelse(all(what == c("projection", "alpha", "H0", "all")) |
           all(what != c("projection", "alpha", "H0", "all")),
         what <- "all",
         what <- what)
  #Choose "parametric" if all or no simulation methods for alpha are choosen
  ifelse(all(asim == c("ordinary", "parametric")) |
           all(asim != c("ordinary", "parametric")),
         ret$asim <- "parametric",
         ret$asim <- asim)
  #Check number of bootstraps
  if(nboot < 1)
    stop("Increase 'nboot', the numbers of bootstrap replicates desired")
  #Create bootno vector
  bootno <- 1 : nboot
  #Keep the number of bootstraps
  ret$nboot <-  nboot
  #Check that 'H0exp' is a list with two elements
  if(!is.list(H0exp) | length(H0exp) != 2)
    stop("'H0exp' should be a list of two elements")
  #Check H0exp and split the expected values into a and M components
  #a
  aexp <- H0exp[[1]]
  if(!(is.vector(H0exp[[1]]) & length(H0exp[[1]]) == x$npar) & !is.null(H0exp[[1]]))
    stop(sprintf("First element of 'H0exp' should either be a vector of length %s or 'NULL'",
                 sQuote(x$npar)))
  #M
  Mexp <- H0exp[[2]]
  if(!is.null(H0exp[[2]]))
    if(is.matrix(H0exp[[2]]))
    {
      if(any(dim(H0exp[[2]]) != c(x$npar, x$npar)))
        stop(sprintf("Second element of 'H0exp' should either be a matrix with dimension '(%s)' or 'NULL'",
                     paste(c(x$npar, x$npar), collapse = ",")))
    }
  else
  {
    stop(sprintf("Second element of 'H0exp' should either be a matrix with dimension '(%s)' or 'NULL'",
                 paste(c(x$npar, x$npar), collapse = ",")))
  }
  #Create the element which takes the optim running time
  if(what == "alpha" | what == "all" | (what == "H0" & any(H0con %in% "ds") & !is.null(Mexp)))
    ret$optim.time <- NA
  #Collect some variables from the model
  ret$uage <- x$uage
  ret$nage <- x$nage
  ret$npar <- x$npar
  ret$nyear <- x$nyear
  ret$uyear <- x$uyear
  #Bootstrap replicates of l (the projection matrix)
  if(what == "projection" | what == "all")
  {
    #Create a matrix with randomly selected individuals (their row number) as
    #rows. Resampling with replacement.
    #If nage > 1. Age ratios are preserved
    if(x$nage > 1)
    {
      samp <- do.call(cbind, sapply(x$uage, function(a) matrix(
        sample(x = c(1 : x$nobs)[x$data.set$age == a], size = x$nobs.age[
          x$uage == a] * nboot, replace = TRUE), ncol = x$nobs.age[x$uage == a])))
    }
    else
    {
      samp <- matrix(sample(x = 1 : x$nobs, size = x$nobs * nboot,
                            replace = TRUE), ncol = x$nobs)
    }
    #l (projection matrix) resampling
    #Collect the estimate from the model
    ret$l <- x$l
    #Resample
    ret$lboot <- lapply(bootno, function(i) promat(procomp(x$data.set[samp[i, ], ],
                                                           uage = x$uage), nage = x$nage))
    #lambda, u (stable age distribution) and v (reproductive values) resampling
    #Collect the estimates from the model
    ret$lambda <- x$lambda
    ret$u <- x$u
    ret$v <- x$v
    #Resample
    ret <- c(ret, list(luvboot = do.call(rbind, lapply(bootno, function(i) unlist(
      eigenl(ret$lboot[[i]]))))))
  }
  #Bootstrap replicates of the alpha parameters and associated parameters
  if(what == "alpha" | what == "all" | (what == "H0" & sig.dj == TRUE) & any(H0con != "nds"))
  {
    #sigma2.dj resampling
    #Collect the estimate from the model
    ret$sigma2.dj <- x$sigma2.dj
    #Resample (get a matrix with npar columns, one for each age)
#    ret$djboot <- mapply(function(a, b) {rnorm(n = nboot, mean = a, sd = b)},
#                         x$sigma2.dj, x$sigma2.dj.sd)
    ret$djboot <- mapply(function(a, b) {rgamma(n = nboot, shape = (a^2)/(b^2), rate = a/(b^2))},
                         x$sigma2.dj, x$sigma2.dj.sd)
    #Assign name to djboot
    ifelse(nboot > 1,
           colnames(ret$djboot) <- paste("sigma2.d", x$uage, sep = ""),
           names(ret$djboot) <- paste("sigma2.d", x$uage, sep = ""))
    #If sigma2.dj contains negative values replace these by zero
#    for(i in 1 : x$nage)
#      ret$djboot[ret$djboot[, i] < 0, ] <- 0
    #sigma.d resampling
    #Collect the estimate from the model
    ret$sigma2.d <- x$sigma2.d
    #Resample
    ret$dboot <- c(ret$djboot %*% x$u)
    #Include bootstraps of sigma2.dj in the bootstraps of alpha and M
    if(sig.dj == TRUE)
    {
      #Calculate the scaling-and-weighting (snw) coefficient (combining resampled
      #age specific demographic variances (sigma2.dj) and stable age distribution(u))
      snw <- t(t(ret$djboot) * x$u^2)
      #At "resampling". Scale Ajt.us by a new bootstrap replicate
      #of sigma2.dj and weight by u^2 for each resampling of sigma2.dj. The
      #last two steps are combined in the "snw" varible
      ifelse(nboot > 1,
             ret$Atboot <- lapply(bootno, function(i)
             {
               lapply(1 : x$nyear, function(j)
               {
                 Reduce('+', mapply('*', sapply(x$Ajt.us, '[', j),
                                    as.list(snw[i, ]), SIMPLIFY = FALSE))
               }
               )
             }
             ),
             ret$Atboot <- lapply(bootno, function(i)
             {
               lapply(1 : x$nyear, function(j)
               {
                 Reduce('+', mapply('*', sapply(x$Ajt.us, '[', j),
                                    as.list(snw[i]), SIMPLIFY = FALSE))
               }
               )
             }
             )
      )
    }
    #Exclude bootstraps of sigma2.dj from the bootstraps of alpha and M
    else
    {
      #At "resampling" without including bootstrapped sigma2.dj. I.e.
      #there is no bootstrap, just set up the estimated At like Atboot.
      ret$Atboot <- rep(list(x$At), nboot)
    }
    #at resampling
    if(ret$asim == "ordinary")
    {
      #ORDINARY
      #Create a matrix with randomly selected years (their row number) as
      #rows. Resampling with replacement
      samp <- matrix(sample(x = 1 : x$nyear, size = x$nyear * nboot, replace = TRUE),
                     ncol = x$nyear)
      #at resampling
      ret$atboot <- lapply(bootno, function(i)'['(x$at, samp[i, ]))
    }
    else
    {
      #PARAMETRIC
      #at resampling from the multinormal distribution. Simulate each at
      #vector one year at the time for each of n bootstraps. For a given
      #bootstrap a given year use At + M as the covariance matrix. Thus,
      #this is specific for each year, but fluctuates within a year
      #between bootstraps if uncertainty in sigma2.dj is included in the
      #estimate of At (M is our estimate and is constrant). The mean is
      #our estimated temporal alpha estimate and is constant.
      ret$atboot <- lapply(bootno, function(i)
      {
        lapply(1 : x$nyear, function(j)
        {
          rmnorm(n = 1, mean = c(x$aM), varcov = (ret$Atboot[[i]][[j]] + x$M))
        }
        )
      }
      )
    }
    #Do not bootstrap M and alpha if only a bootstrap under the null
    #hypothesis of no fluctuating selection with uncertainty in the
    #demographic variance is desired
    if(what == "alpha" | what == "all" | (what == "H0" & "fs" %in% H0con))
    {
      #M and aM resampling. Bootstrapping of fluctuating selection
      #Start timer
      ret$optim.time <- proc.time()[3]
      #Collect the estimate from the model
      ret$aM <- x$aM
      ret$M <- x$M
      #Resampling
      #Estimate fluctuating selection for the first bootstrap replicate
      ret$Mboot <- vector("list", nboot)
      ret$aMboot <- vector("list", nboot)
      for (i in bootno)
      {
        fs <- fs(At = ret$Atboot[[i]], at = ret$atboot[[i]], npar = x$npar,
                 nyear = x$nyear, method = method, control = control, ...)
        ret$Mboot[[i]] <- fs$M
        ret$aMboot[[i]] <- fs$aM
      }
      #Stop timer
      ret$optim.time <- c(proc.time()[3] - ret$optim.time)
      #E(at|?t) resampling
      ret$atCboot <- mapply(function(a, b, c, d)
      {atCfn(aM = a, M = b, At = c, at = d)},
                            ret$aMboot, ret$Mboot, ret$Atboot, ret$atboot, SIMPLIFY = FALSE)
    }
    #a|M=0 and A|M=0 resampling. Bootstrap alphas and
    #alphas covariance matrix under the assumption of
    #no fluctuating selection (M = 0)
    #Collect the estimates from the model
    ret$Anf <- x$Anf
    ret$anf <- x$anf
    #Resampling
    ret$Anfboot <- vector("list", nboot)
    ret$anfboot <- vector("list", nboot)
    for (i in bootno)
    {
      nfs <- nfs(At = ret$Atboot[[i]], at = ret$atboot[[i]], npar = x$npar,
                 nyear = x$nyear)
      ret$Anfboot[[i]] <- nfs$Anf
      ret$anfboot[[i]] <- nfs$anf
    }
    #sigma2.e resampling normal and with correction for sampling error
    #Collect the estimates from the model
    ret$sigma2.e <- x$sigma2.e
    #Normal
    ret$eboot <- mapply(function(a)var(sapply(a, '[', 1)), ret$atboot)
  }
  #Bootstrap replicates of the alpha vector and the M - matrix under the
  #specified null hypothesis
  if(what == "H0" | what == "all" & (!is.null(aexp) | !is.null(Mexp)))
  {
    #1. H0: a=0|M
    if(!is.null(aexp) & "fs" %in% H0con)
    {
      #Check if aM and M estimates have been collected earlier
      if(is.null(ret$aM))
      {
        #Collect the estimate from the model
        ret$aM <- x$aM
        ret$M <- x$M
      }
      #Resample alpha vectors under the null hypothesis
      #With demographic variance
      if(sig.dj == TRUE)
      {
        #Set up the matrix with resampled alphas
        ret$H0aMboot <- matrix(rep(NA, x$npar * nboot), ncol = x$npar)
        #Insert names
        colnames(ret$H0aMboot) <- colnames(x$aM)
        #Bootstrap, choose a new M for each resampling to include uncertainty
        #in sigma2.dj which was resampled earlier
        for(i in bootno)
        {
          # alphas are assumed multinormally distributed with covariance matrix M,
          # see Appendix A.
          ret$H0aMboot[i, ] <- rmnorm(n = 1, mean = aexp, varcov = ret$Mboot[[i]])
        }
      }
      #Without demographic variance
      else
      {
        #Bootstrap without including uncertainty in sigma2.dj
        ret$H0aMboot <- rmnorm(n = nboot, mean = aexp, varcov = x$M)
        #Insert names
        colnames(ret$H0aMboot) <- colnames(x$aM)
      }
    }
    #2. H0: a=0|M=0
    if(!is.null(aexp) & "nfs" %in% H0con)
    {
      #Check if aM, M, anf and Anf estimates have been collected earlier
      if(is.null(ret$aM))
      {
        #Collect the estimate from the model
        ret$aM <- x$aM
        ret$M <- x$M
        ret$anf <- x$anf
        ret$Anf <- x$Anf
      }
      #Resample a vectors under the null hypothesis
      #With demographic variance
      if(sig.dj == TRUE)
      {
        #Set up the matrix with resampled alphas
        ret$H0anfboot <- matrix(rep(NA, x$npar * nboot), ncol = x$npar)
        #Insert names
        colnames(ret$H0anfboot) <- colnames(x$aM)
        #Bootstrap, choose a new Anf for each resampling to include uncertainty
        #in sigma2.dj which was resampled earlier
        for(i in bootno)
        {
          ret$H0anfboot[i, ] <- rmnorm(n = 1, mean = aexp, varcov = ret$Anfboot[[i]])
        }
      }
      #Without demographic variance
      else
      {
        #Bootstrap without including uncertainty in sigma2.dj
        ret$H0anfboot <- rmnorm(n = nboot, mean = aexp, varcov = x$Anf)
        #Insert names
        colnames(ret$H0anfboot) <- colnames(x$aM)
      }
    }
    #3. H0: M=0|a
    if(!is.null(Mexp) & "ds" %in% H0con)
    {
      #Check if aM, M, anf and Anf estimates have been collected earlier
      if(is.null(ret$aM))
      {
        #Collect the estimate from the model
        ret$aM <- x$aM
        ret$M <- x$M
        ret$anf <- x$anf
        ret$Anf <- x$Anf
      }
      #Resample at vectors from the null hypothesis distribution
      #With demographic variance
      if(sig.dj == TRUE)
      {
        ret$H0atnfboot <- lapply(bootno, function(i)
        {
          lapply(1 : x$nyear, function(j)
          {
            rmnorm(n = 1, mean = c(x$anf), varcov = ret$Atboot[[i]][[j]])
          }
          )
        }
        )
      }
      #Without demographic variance
      else
      {
        ret$H0atnfboot <- lapply(bootno, function(i)
        {
          lapply(1 : x$nyear, function(j)
          {
            rmnorm(n = 1, mean = c(x$anf), varcov = x$At[[j]])
          }
          )
        }
        )
      }
      #For each of t resampled at vectors calculate M under the null hypothesis
      #Start timer
      H0M.ot <- proc.time()[3]
      #Set up the list which takes the bootstraps
      ret$H0Mnfboot <- vector("list", nboot)
      for (i in bootno)
      {
        fs <- fs(At = x$At, at = ret$H0atnfboot[[i]], npar = x$npar,
                 nyear = x$nyear, method = method, control = control, ...)
        ret$H0Mnfboot[[i]] <- fs$M
      }
      #Stop timer
      H0M.ot <- c(proc.time()[3] - H0M.ot)
      if(is.na(ret$optim.time))
      {
        ret$optim.time <- H0M.ot
      }
      else
      {
        ret$optim.time <- ret$optim.time + H0M.ot
      }
    }
    #4. H0: M=0|a=0
    if(!is.null(Mexp) & "nds" %in% H0con)
    {
      #No 4. might inflate the risk of Type I error if the assumption is
      #wrong, i.e. a != 0, but there is not enough power to detect it.
      warning("Simulations under the null hypothesis M=0|a=0 inflates risk of Type I error if the assumption a = 0 is wrong, thus it has not been implemented")
    }
  }
  #Calculate total running time
  ret$running.time <- proc.time()[3] - ret$running.time
  #Output
  class(ret) <- "boot.lmf"
  ret
}
