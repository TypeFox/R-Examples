summary.boot.lmf <-
function(object,
                             ret.bootstraps = FALSE,
                             ...)
{
  #Rename the object
  z <- object
  #Collect call
  ans <- z[c("call")]
  #Extract number of bootstraps
  ans$nboot <- z$nboot
  #Create summary of projection matrix and associated components
  if(!is.null(z$lboot))
  {
    #Projection matrix
    #Estimate
    ans$lest <- z$l
    #Set up as array
    lt <- array(Reduce(cbind, z$lboot),
                dim = c(dim(z$l)[1], dim(z$l)[2], z$nboot))
    #Mean from bootstrap
    ans$lboot.mean <- apply(lt, 1:2, mean, na.rm = TRUE)
    #Bias, i.e. estimate - mean from bootstrap
    ans$lbias <- ans$lest - ans$lboot.mean
    #Standard deviation from bootstrap
    ans$lboot.sd <- apply(lt, 1:2, sd, na.rm = TRUE)
    #Lambda, u and v
    #Estimates
    luvest <- c(z$lambda, z$u, z$v)
    #Means from bootstraps
    luvboot.mean <- colMeans(z$luvboot, na.rm = TRUE)
    #Bias, i.e. estimates - means from bootstraps
    luvbias <- luvboot.mean - luvest
    #Standard deviations from bootstraps
    luvboot.sd <- apply(z$luvboot, 2, sd, na.rm = TRUE)
    #Put together as a data frame
    ans$luv <- cbind(estimate = luvest, boot.mean = luvboot.mean,
                     bias = luvbias, boot.sd = luvboot.sd)
  }
  #Create summary statistics of the remaining parameters (all related
  #to the alpha estimates)
  if(!is.null(z$aMboot))
  {
    #Variance components
    #Environmental
    #Normal
    #Estimate
    eest <- z$sigma2.e
    #Mean of bootstraps
    eboot.mean <- mean(z$eboot)
    #Bias, i.e. estimate - mean from bootstraps
    ebias <- eboot.mean - eest
    #Standard deviation from bootstraps
    eboot.sd = sd(z$eboot)
    #Put together as a data frame
    ans$sigma2.e <- data.frame(estimate = eest, boot.mean = eboot.mean,
                               bias = ebias, boot.sd = eboot.sd)
    #Demographic
    #Estimate
    djest <- unlist(z$sigma2.dj)
    dest <- z$sigma2.d
    #Mean of bootstraps
    djboot.mean <- as.vector(apply(z$djboot, 2, mean))
    dboot.mean <- mean(z$dboot)
    #Bias, i.e. estimate - mean from bootstraps
    djbias <- djboot.mean - djest
    dbias <- dboot.mean - dest
    #Standard deviation from bootstraps
    djboot.sd <- as.vector(apply(z$djboot, 2, sd))
    dboot.sd <- sd(z$dboot)
    #Put together as a data frame
    ans$sigma2.dd <- data.frame(age = c(z$uage, "(total)"),
                                estimate = c(djest, dest), boot.mean = c(djboot.mean, dboot.mean),
                                bias = c(djbias, dbias), boot.sd = c(djboot.sd, dboot.sd))
    #alphas (aM) - Fluctuating selection
    #Estimate
    aMest <- c(z$aM)
    #Set up as data frame with alpha estimates in separate columns
    aMt <- do.call(rbind, z$aMboot)
    #Mean of bootstraps
    aMboot.mean <- colMeans(aMt, na.rm = TRUE)
    #Bias, i.e. estimate - mean from bootstraps
    aMbias <- aMboot.mean - aMest
    #Standard deviation from bootstraps
    aMboot.sd = apply(aMt, 2, sd, na.rm = TRUE)
    #Put together as a data frame
    ans$aM <- cbind(estimate = aMest, boot.mean = aMboot.mean,
                    bias = aMbias, boot.sd = aMboot.sd)
    #alphas covariance matrix (M) - Fluctuating selection
    #Estimate
    ans$Mest <- z$M
    #Set up as array
    Mt <- array(Reduce(cbind, z$Mboot),
                dim = c(z$npar, z$npar, z$nboot))
    #Mean from bootstrap
    ans$Mboot.mean <- apply(Mt, 1:2, mean, na.rm = TRUE)
    dimnames(ans$Mboot.mean) <- dimnames(ans$Mest)
    #Bias, i.e. estimate - mean from bootstrap
    ans$Mbias <- ans$Mest - ans$Mboot.mean
    #Standard deviation from bootstrap
    ans$Mboot.sd <- apply(Mt, 1:2, sd, na.rm = TRUE)
    dimnames(ans$Mboot.sd) <- dimnames(ans$Mest)
    #alphas (a(M=0)) - No fluctuating selection
    #Estimate
    anfest <- c(z$anf)
    #Set up as data frame with alpha estimates in separate columns
    anft <- do.call(rbind, z$anfboot)
    #Mean of bootstraps
    anfboot.mean <- colMeans(anft)
    #Bias, i.e. estimate - mean from bootstraps
    anfbias <- anfboot.mean - anfest
    #Standard deviation from bootstraps
    anfboot.sd = apply(anft, 2, sd)
    #Put together as a data frame
    ans$anf <- cbind(estimate = anfest, boot.mean = anfboot.mean,
                     bias = anfbias, boot.sd = anfboot.sd)
    #alphas covariance matrix (A) - No fluctuating selection
    #Estimate
    ans$Anfest <- z$Anf
    #Set up as array
    Anft <- array(Reduce(cbind, z$Anfboot),
                  dim = c(z$npar, z$npar, z$nboot))
    #Mean from bootstrap
    ans$Anfboot.mean <- apply(Anft, 1:2, mean)
    dimnames(ans$Anfboot.mean) <- dimnames(ans$Anfest)
    #Bias, i.e. estimate - mean from bootstrap
    ans$Anfbias <- ans$Anfest - ans$Anfboot.mean
    #Standard deviation from bootstrap
    ans$Anfboot.sd <- apply(Anft, 1:2, sd)
    dimnames(ans$Anfboot.sd) <- dimnames(ans$Anfest)
  }
  #Create matrix with the tests of significance for the hypotheses which
  #have been runned.
  #H0: a=0|M
  if(!is.null(z$H0aMboot))
  {
    #Extract estimates with se, test-statistic and p-value
    est.aM <- as.vector(z$aM)
    boot.se.aM <-  apply(z$H0aMboot, 2, sd, na.rm = TRUE)
    n <- colSums(t(t(z$H0aMboot) > as.vector(est.aM)), na.rm = TRUE)
    nsucc.aM <- sapply(n, function(a) min(a, (sum(!is.na(z$H0aMboot[, 1])) - a)))
    p <- n / sum(!is.na(z$H0aMboot[, 1]))
    pval.aM <- sapply(p, function(a) min((2 * a), 2 * (1 - a)))
    ans$coefficients.H0aMboot <- cbind(est.aM, boot.se.aM, nsucc.aM, pval.aM)
    dimnames(ans$coefficients.H0aMboot) <- list(dimnames(z$aM)[[2]],
                                                c("Estimate", "Std. Error", "# successes", "Pr(>|(a(M))|)"))
    if(sum(!is.na(z$H0aMboot[, 1])) < z$nboot){
      H0aMboot.nboot <- sum(!is.na(z$H0aMboot[, 1]))
      msg <- sprintf("Only %s bootstraps of temporal mean coefficients under the null hypothesis of fluctuating selection obtained due to lack of convergence. P-values corrected", H0aMboot.nboot)
      warning(msg, domain = NA)
    }
      
  }
  #H0: a=0|M=0
  if(!is.null(z$H0anfboot))
  {
    #Extract estimates with se, test-statistic and p-value
    est.anf <- as.vector(z$anf)
    boot.se.anf <-  apply(z$H0anfboot, 2, sd)
    n <- colSums(t(t(z$H0anfboot) > as.vector(est.anf)))
    nsucc.anf <- sapply(n, function(a) min(a, (z$nboot - a)))
    p <- n / z$nboot
    pval.anf <- sapply(p, function(a) min((2 * a), 2 * (1 - a)))
    ans$coefficients.H0anfboot <- cbind(est.anf, boot.se.anf, nsucc.anf, pval.anf)
    dimnames(ans$coefficients.H0anfboot) <- list(dimnames(z$anf)[[2]],
                                                 c("Estimate", "Std. Error", "# successes", "Pr(>|(a(M=0))|)"))
  }
  #H0: M=0|a
  if(!is.null(z$H0atnfboot))
  {
    #Extract names for the M matrix
    Mnames <- vector("list", z$npar)
    #The for loop create names which are of the type (Intercept)-(Intercept)
    #for only the six unique pairwise combinations
    for(i in 1 : z$npar)
      Mnames[[i]] <- paste(colnames(z$M)[i], colnames(z$M)[i : z$npar], sep = "-")
    #Extract estimates with se, test-statistic and p-value
    est.M <- z$M[lower.tri(z$M, diag = TRUE)]
    #Insert M|a(M=0) bootstraps
    H0Mnfboot <- t(sapply(z$H0Mnfboot, function(a){ a[lower.tri(a, diag = TRUE)]}))
    #Insert names
    colnames(H0Mnfboot) <- paste(unlist(Mnames), "(M|a(M=0))", sep = "")
    #Calculate standard deviation
    boot.se.M <-  apply(H0Mnfboot, 2, sd, na.rm = TRUE)
    #Calculate the number of simulated numbers above the estimated ones
    n <- colSums(t(t(H0Mnfboot) > as.vector(est.M)), na.rm = TRUE)
    #Calculate the number of successes (numbers of times the estimate
    #is larger than the simulated expected values under the null hypothesis)
    nsucc.M <- sapply(n, function(a) min(a, (sum(!is.na(H0Mnfboot[, 1])) - a)))
    #Calculate the proportion n to number of bootstraps
    p <- n / sum(!is.na(H0Mnfboot[, 1]))
    #Calculate the p-value
    pval.M <- sapply(p, function(a) min((2 * a), 2 * (1 - a)))
    #Set up matrix
    ans$coefficients.H0Mnfboot <- cbind(est.M, boot.se.M, nsucc.M, pval.M)
    #Give matrix names
    dimnames(ans$coefficients.H0Mnfboot) <- list(unlist(Mnames),
                                                 c("Estimate", "Std. Error", "# successes", "Pr(>|(M|a)|)"))
    if(sum(!is.na(H0Mnfboot[, 1])) < z$nboot){
      H0Mnfboot.nboot <- sum(!is.na(H0Mnfboot[, 1]))
      msg <- sprintf("Only %s bootstraps of temporal covariance matrix under the null hypothesis of directional, but no fluctuating selection obtained due to lack of convergence. P-values corrected", H0Mnfboot.nboot)
      warning(msg, domain = NA)
    }
  }
  #Return bootstraps in a clear way if desired
  if(ret.bootstraps)
  {
    if(!exists(as.character(substitute(Mnames))))
    {
      #Extract names for the M matrix
      Mnames <- vector("list", z$npar)
      #The for loop create names which are of the type (Intercept)-(Intercept)
      #for only the six unique pairwise combinations
      for(i in 1 : z$npar)
        Mnames[[i]] <- paste(colnames(z$M)[i], colnames(z$M)[i : z$npar], sep = "-")
    }
    #Check if projection matrix is bootstrapped and return
    if(!is.null(z$lboot))
    {
      #l, lambda, u and v
      #Set up data frame which takes the non-zero part of the projection
      #matrix components and extract components. First check if
      #individuals in age class c (the final age class) have a
      #survival probability (i.e. survival is not zero)
      if(z$l[z$nage, z$nage] == 0)
      {
        #Set up matrix (when == 0)
        lboot <- matrix(rep(NA, (z$nage * 2 - 1) * z$nboot), ncol = z$nage * 2 - 1)
        #Assign names (when == 0)
        colnames(lboot) <- c(paste("f", 1 : z$nage, sep = ""), paste("s", 1 : (z$nage - 1), sep = ""))
        #Extract survivals (when == 0)
        if(z$nage == 2)
        {
          lboot[, z$nage + 1] <- sapply(z$lboot, '[', z$nage, 1)
        }
        else
        {
          for(i in 1 : z$nage - 1)
            lboot[, z$nage + i] <- sapply(z$lboot, function(a) diag(a[-1, ])[i])
        }
      }
      else
      {
        #Set up matrix (when != 0)
        lboot <- matrix(rep(NA, z$nage * 2 * z$nboot), ncol = z$nage * 2)
        #Assign names (when != 0)
        colnames(lboot) <- c(paste("f", 1 : z$nage, sep = ""), paste("s", 1 : z$nage, sep = ""))
        #Extract survivals (when != 0)
        if(z$nage == 2)
        {
          lboot[, z$nage + 1] <- sapply(z$lboot, '[', z$nage, 1)
          lboot[, z$nage + 2] <- sapply(z$lboot, '[', z$nage, 2)
        }
        else
        {
          for(i in 1 : z$nage - 1)
            lboot[, z$nage + i] <- sapply(z$lboot, function(a) diag(a[-1, ])[i])
          lboot[, z$nage * 2] <- sapply(z$lboot, '[', z$nage, z$nage)
        }
      }
      #Extract fecundities (This is the same procedure whether survival
      #of age class is zero or not)
      for(i in 1 : z$nage)
        lboot[, i] <- sapply(z$lboot, '[', 1, i)
      #Insert the components of the projection matrix into the
      #projection bootstrap data frame in 'ans' list
      ans$lluvboot <- cbind(lboot, z$luvboot)
    }
    #Check if sigma.e and related parameters is bootstrapped and return
    if(!is.null(z$eboot))
    {
      #sigma2.dj, sigma2.d, sigma2.e and sigma2.eC
      #Set up matrix which take the bootstrapped variance components
      ans$deboot <- cbind(z$djboot, z$dboot, z$eboot)
      #Insert names
      colnames(ans$deboot) <- c(colnames(z$djboot), "sigma2.d", "sigma2.e")
    }
    #Check if at and related parameters is bootstrapped and return
    if(!is.null(z$atboot))
    {
      #at and At
      #Set up matrix which take the bootstrapped yearly alphas and
      #covariance components (number of columns is npar for alpha estimates
      #and npar * (npar + 1) / 2 for the covariance matrix and two more
      #for year and bootstrap number columns
      ans$atAboot <- matrix(rep(NA,
                                (z$npar + z$npar * (z$npar + 1) / 2 + 2) * z$nyear * z$nboot),
                            ncol = z$npar + z$npar * (z$npar + 1) / 2 + 2)
      #Insert names
      #Insert names
      colnames(ans$atAboot) <- c("bootno", "year",
                                 paste(colnames(z$aM), "(at)", sep = ""),
                                 paste(unlist(Mnames), "(At)", sep = ""))
      #Insert 'bootno' and 'year'
      ans$atAboot[, 1] <- rep(1 : z$nboot, each = z$nyear)
      ans$atAboot[, 2] <- rep(z$uyear, times = z$nboot)
      #Insert alpha (at) bootstraps
      ans$atAboot[, 3 : (2 + z$npar)] <- do.call(rbind,
                                                 mapply(function(a) do.call(rbind, a), z$atboot, SIMPLIFY = FALSE))
      #Insert At bootstraps
      ans$atAboot[, (3 + z$npar) :
                    ((2 + z$npar) + z$npar * (z$npar + 1) / 2)] <- do.call(rbind,
                                                                           mapply(function(a) t(sapply(a, function(b){
                                                                             b[lower.tri(b, diag = TRUE)]})), z$Atboot, SIMPLIFY = FALSE))
    }
    #Check if aM and related parameters is bootstrapped and return
    if(!is.null(z$aMboot))
    {
      #aM and M
      #Set up matrix which take the bootstrapped temporal alphas and
      #the temporal covariance components
      ans$aMMboot <- matrix(rep(NA,
                                (z$npar + z$npar * (z$npar + 1) / 2) * z$nboot),
                            ncol = z$npar + z$npar * (z$npar + 1) / 2)
      #Insert names
      colnames(ans$aMMboot) <- c(paste(colnames(z$aM), "(a(M))", sep = ""),
                                 paste(unlist(Mnames), "(M)", sep = ""))
      #Insert alpha (aM) bootstraps
      ans$aMMboot[, 1 : z$npar] <- do.call(rbind, z$aMboot)
      #Insert M bootstraps
      ans$aMMboot[, (1 + z$npar) :
                    (z$npar + z$npar * (z$npar + 1) / 2)] <-
        t(sapply(z$Mboot, function(a){ a[lower.tri(a, diag = TRUE)]}))
      #atC
      #Set up matrix which take the bootstrapped yearly corrected alphas
      ans$atCboot <- matrix(rep(NA, (z$npar + 2) * z$nyear * z$nboot),
                            ncol = z$npar + 2)
      #Insert names
      colnames(ans$atCboot) <- c("bootno", "year",
                                 paste(colnames(z$aM), "(atC)", sep = ""))
      #Insert 'bootno' and 'year'
      ans$atCboot[, 1] <- rep(1 : z$nboot, each = z$nyear)
      ans$atCboot[, 2] <- rep(z$uyear, times = z$nboot)
      #Insert alpha (atC) bootstraps
      ans$atCboot[, 3 : (2 + z$npar)] <- do.call(rbind,
                                                 mapply(function(a) do.call(rbind, a), z$atCboot, SIMPLIFY = FALSE))
    }
    #Check if anf and related parameters have been bootstrapped and return
    if(!is.null(z$anfboot))
    {
      #anf and Anf
      #Set up matrix which take the bootstrapped temporal alphas and
      #the temporal covariance components assuming no fluctuating
      #selection
      ans$anfAboot <- matrix(rep(NA,
                                 (z$npar + z$npar * (z$npar + 1) / 2) * z$nboot),
                             ncol = z$npar + z$npar * (z$npar + 1) / 2)
      #Insert names
      colnames(ans$anfAboot) <-
        c(paste(colnames(z$aM), "(a(M=0))", sep = ""),
          paste(unlist(Mnames), "(A)", sep = ""))
      #Insert alpha (anf) bootstraps
      ans$anfAboot[, 1 : z$npar] <- do.call(rbind, z$anfboot)
      #Insert Anf bootstraps
      ans$anfAboot[, (1 + z$npar) :
                     (z$npar + z$npar * (z$npar + 1) / 2)] <-
        t(sapply(z$Anfboot, function(a){ a[lower.tri(a, diag = TRUE)]}))
    }
    #Check if bootstraps under a null hypothesis is found and return
    #H0: a = 0 | M
    if(!is.null(z$H0aMboot))
    {
      #Keep the bootstrapped alpha estimates under the nullhypothesis
      #given that there is fluctuating selection
      ans$H0aMboot <- z$H0aMboot
      #Insert names
      colnames(ans$H0aMboot) <- c(paste(colnames(z$aM), "(a=0|M)", sep = ""))
    }
    #H0: a = 0 | M = 0
    if(!is.null(z$H0anfboot))
    {
      #Keep the bootstrapped alpha estimates under the nullhypothesis
      #given that there is no fluctuating selection
      ans$H0anfboot <- z$H0anfboot
      #Insert names
      colnames(ans$H0anfboot) <- c(paste(colnames(z$aM), "(a=0|M=0)", sep = ""))
    }
    #H0: M = 0 | a(M=0)
    if(!is.null(z$H0atnfboot))
    {
      #Set up matrix which take the bootstrapped yearly alphas (at|M = 0)
      #under the nullhypothesis
      ans$H0atnfboot <- matrix(rep(NA, (z$npar + 2) * z$nyear * z$nboot),
                               ncol = z$npar + 2)
      #Insert names
      colnames(ans$H0atnfboot) <- c("bootno", "year",
                                    paste(colnames(z$aM), "(at|a)", sep = ""))
      #Insert 'bootno' and 'year'
      ans$H0atnfboot[, 1] <- rep(1 : z$nboot, each = z$nyear)
      ans$H0atnfboot[, 2] <- rep(z$uyear, times = z$nboot)
      #Insert alpha (at|M = 0) bootstraps
      ans$H0atnfboot[, 3 : (2 + z$npar)] <- do.call(rbind,
                                                    mapply(function(a) do.call(rbind, a), z$H0atnfboot, SIMPLIFY = FALSE))
      #Set up matrix which take the bootstrapped temporal covariance
      #components (M) bootstrapped under the nullhypothesis given that
      #there is directional selection
      #Insert M|M = 0 bootstraps
      ans$H0Mnfboot <- H0Mnfboot
    }
  }
  #Define class
  class(ans) <- "summary.boot.lmf"
  #Return
  ans
}
