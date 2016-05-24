#
# pcse.R
# function to estimate panel-corrected standard errors
#

pcse <- function(object, groupN, groupT, pairwise=FALSE){
  # 
  # Purpose:
  # To estimate panel-corrected standard errors.
  # 
  # Input arguments:
  # object   = an lm object
  # groupN   = a vector containing the cross-sectional group identifier
  #            for each observation
  # groupT   = a vector containing the time-series identifier for each
  #            observation
  # pairwise = An optional logical flag indicating whether the X's used to
  #            estimate the "middle" matrix should be chosen in a pairwise
  #            fashion or casewise fashion. If pairwise, the correlation between
  #            observations $i$ and $j$ is based on the time periods common to
  #            $i$ and $j$. If casewise, the correlation between observations i
  #            and j is based on the largest rectangular subset of the data,
  #            i.e., $T_i$ = $T_j$ = $T^*$ for all $i$ and $j$ if casewise is
  #            selected.
  #
  # Output values:
  # vcov     = the panel-corrected variance covariance matrix
  # pcse     = panel-corrected standard errors
  # b        = beta values as passed by the 'lm' object
  # tstats   = t-statistics 
  # df       = degrees of freedom
  # pval     = p-values
  # pairwise = logical flag corresponding to pairwise argument input
  # nobs     = number of observations
  # nmiss    = number of missing observations relative to the full balanced panel
  # call     = the call to pcse
  #
  # Exits:
  # pcse will fail with a warning if:
  # 1. the formula object is not of the class 'lm'
  # 2. the time-series and cross-section identifier variables are of differing
  #    length
  # 3. The data is of differing length than the time-series and cross-section
  #    identifier variables
  # 4. There are missing values in the cross-section identifier variable
  # 5. There are missing values in the time-series identifier variable
  # 6. There are more than nCS*nTS rows in the data.
  # 7. If casewise selection is used and all observations are deleted.
  # 8. If pairwise selection is used and there is a cross-section unit with
  #    no observations in common with another cross-section unit.
  # 


  # Get call to pcse.
  mc <- match.call()

  # Input Checks

  # Check the formula object is 'lm.'
  check <- class(object)
  if (!("lm" %in% check)){
    stop("Formula object must be of class 'lm'.")
  }


  # Get identifier variables in form we like --
  # (can have data$groupT, groupT or "groupT") supplied originally
  extractLMvar <- function(lm.result, groupvar){
    # extract named argument
    isNamed <- try(inherits(groupvar, "character"), silent=TRUE)
    if (inherits(isNamed, "try-error")){
      # handle unquoted name
      groupvar <- as.character(match.call()$groupvar)
      isNamed <- TRUE
    }
    if (isNamed & length(groupvar) == 1){
      groupvar <- model.frame(lm.result)[[groupvar]]
    }
    groupvar
  }
  groupT <- extractLMvar(object, groupT)
  groupN <- extractLMvar(object, groupN)
  # check that groupT and groupN were found in the data, global namespace, or model
  check <- !is.null(groupT)
  if (!check){
      stop(paste(deparse(mc$groupT), "is not found."))
  }
  check <- !is.null(groupN)
  if (!check){
      stop(paste(deparse(mc$groupN), "is not found."))
  }
      
  
  # Check that the time-series and cross-section identifier variables
  # are the same length.
  check <- length(groupN) == length(groupT)
  if (!check){
    stop("Length of groupT and groupN be of equal length.")
  }

  # Check that time-series and cross-section identifier variables
  # are the same length as the using data.
  check <- length(groupN) == dim(model.matrix(object))[1]
  if (!check){
    stop("Length of groupN and groupT must equal nrows of using data.")
  }

  # Check that there are not any missing values in the cross-section
  # identifier variable.
  check <- is.na(groupN)
  if (any(check)){
    stop("There must not be any missing values in the CS groupN!")
  }

  # Check that there are not any missing values in the time-series
  # identifier variable.
  check <- is.na(groupT)
  if (any(check)){
    stop("There must not be any missing values in the TS groupT!")
  }


  # Check that there are not more than nCS*nTS rows in the data.
  nCS     <- length(na.omit(unique(groupN))) # number of cross-sectional units
  nTS     <- length(na.omit(unique(groupT))) # number of time-series units
  check <- nCS*nTS >= dim(model.matrix(object))[1]
  if (!check){
    stop("There cannot be more than nCS*nTS rows in the using data!")
  }

  # Make factors numeric for easier manipulation.
  if ("factor" %in% class(groupN)){
    groupN <- as.numeric(groupN)
  }
  if ("factor" %in% class(groupT)){
    groupT <- as.numeric(groupT)
  }
 
  # Check for balanced data.
  units <- unique(groupN) 
  units <- na.omit(units) # Unique cross-sectional units.
  nCS   <- length(units)  # Number of unique cross-sectional units.
  time  <- unique(groupT)
  time  <- na.omit(time)  # Unique time-series units.
  nTS   <- length(time)   # Number of unique time-series units.

  check <- 0
  for (i in 1:nCS){       # Loop over nCS and check if number of observations
                          # per panel are equal across panels.
    if (sum(groupN  == units[i]) == nTS){
      check <- check + 1
    }
    else{
      check <- check
    }
  }
  flag <- ifelse(check == nCS, TRUE, FALSE) # flag for balanced data.
  # END CHECKS

  # Data manipulation.

  # Make data frame with [CS Id, TS Id, Residuals, Model matrix from 'lm' call]
  using       <- data.frame(groupN = groupN, groupT = groupT)
  using       <- na.omit(using) 
  using$resid <- resid(object)  
  using       <- data.frame(using, model.matrix(object)) 
  

  # Sort data by cross-sectional ID then time-series ID
  ord <- order(using$groupN, using$groupT)
  using <- using[ord, ]

  # Get list of unique cross-sectional and time-series unit identifiers and
  # number of unique units in each.
  units <- unique(using$groupN) # Unique cross-sectional units.
  nCS   <- length(units)        # Number of unique cross-sectional units.
  time  <- unique(using$groupT) # Unique time-series units.
  nTS   <- length(time)         # Number of unique time-series units.

  # calculate avg number of obs per panel.
  avgN <- dim(using)[1] / nCS

  # get largest balanced subset of data.
  brows <- c()          # Place holder for balanced rows.
  for (i in 1:nTS){     # Loop over number of time-series units.
    # Get row numbers corresponding to time-series unit[i] and check that
    # there are nCS rows for time-series unit[i].
    # If so, keep those rows. If not, pass. 
    br    <- which(using$groupT == time[i])
    check <- length(br) == nCS
    if (check){
      brows <- c(brows, br)
    }
  }

  # Pull out the rows that are balanced subsets and re-sort by CS Id then TS Id.
  balanced <- using[brows, ]
  ord      <- order(balanced$groupN, balanced$groupT)
  balanced <- balanced[ord, ]

  # Get list of unique cross-sectional and time-series unit identifiers within
  # the balanced subset of the data and the number of units in each.
  Bunits <- unique(balanced$groupN)
  BnCS   <- length(Bunits)
  Btime  <- unique(balanced$groupT)
  BnTS   <- length(Btime)
 
  # Get rectangular data.
  rect <- using

  # Fill in missing Time and Cross-Section Units.
  rect$groupN <- as.numeric(rect$groupN) 
  Runits      <- unique(rect$groupN)
  Runits      <- na.omit(Runits) # Unique cross-sectional units in the
                                # rectangular data.
 
  # What are we missing in each panel?
  # Loop over unique cross-sectional units and if the panel is not balanced put
  # the CS unit in the first column of the matrix and the number of missing obs
  # in the panel in the second column.
  missN <- matrix(NA, nCS, 2)           # Matrix to hold missing rows #'s
  for (i in 1:nCS){                   
    if (sum(rect$groupN == Runits[i]) != nTS){
      missN[i, 1] <- Runits[i]
      missN[i, 2] <- nTS - (sum(rect$groupN == Runits[i]))
    }
  }

  # Throw out any empty rows (correspond to balanced panel units)
  missN <- na.omit(missN)
  
  missT <- c() # Holder for the missing time units
  tmp   <- c() # temporary holder.

  # If the missing row matrix is not empty, then
  # loop over the number of rows in the missing matrix
  # and pull out the time periods that are missing for that CS unit.
  if (dim(missN)[1] != 0 & dim(missN)[2] != 0){
    for (i in 1:dim(missN)[1]){
      # Time periods that appear in the data for panel[i]
      tt <- time %in% rect$groupT[rect$groupN == missN[i, 1]]
      # Grab the missing time periods
      missT <- c(missT, time[!tt])
      # temporary vector which expands CS unit [i] the number of
      # necessary times 
      tmp <- c(tmp, rep(missN[i, 1], missN[i, 2]))
    }
    missN <- tmp # After looping over all rows, replace with temporary frame
    nM    <- length(missN) # Number of observations to replace

    # Fill in using data with NAs to get N*T rows.
    R <- dim(rect)[1]                # nrows in rectangular data
    C <- dim(rect)[2]                # ncols in rectangular data
    if (R != nCS * nTS){               # if nrows != nCS*nTS
      for (i in (R+1):(nCS * nTS)){        # loop over the extra rows
        rect[i, ] <- rep(NA, C)      # and fill in NAs for those rows
      }
    }
    rect[c((R+1):(R+nM)), 1] <- missN # fill in 1st col of missing rows with
                                     # missing cross-sectional units
    rect[c((R+1):(R+nM)), 2] <- missT # fill in 2nd col of missing rows with
                                     # missing time-series units
  }

  # Sort rectangular data by CS Id then TS Id.
  ord    <- order(rect$groupN, rect$groupT)
  rect   <- rect[ord, ]

  # Get list of unique cross-sectional and time-series unit identifiers within
  # the rectangular data and the number of units in each.
  Runits <- unique(rect$groupN) 
  Runits <- na.omit(Runits)
  RnCS   <- length(Runits)
  Rtime  <- unique(rect$groupT)
  Rtime  <- na.omit(Rtime)
  RnTS   <- length(Rtime)
  # End data manipulation.
  
  # ESTIMATION

  # IF Balanced data, then:
  if (flag){
    # estimate Sigma.hat using whole data and get middle matrix (X' omega X).

    e <- using$resid            # residuals
    # Reshape the vector of residuals into the nCS x nTS matrix
    E <- matrix(e, nCS, nTS, byrow=TRUE)
    # Actually want it to be nTS x nCS.
    E <- t(E)

    # Sigma.hat = (E'E)/nTS
    Sigma.hat <- crossprod(E)/nTS
    # X = data without first 3 cols (CS Id, TS Id, residuals)
    X         <- as.matrix(using[, 4:dim(using)[2]])
    # omega = Sigma.hat kronecker Identity_nT
    omega     <- kronecker(Sigma.hat, diag(1, nTS))
    # middle matrix = X' omega X
    middle    <- t(X) %*% omega %*% X
    # number of observations = number of residuals
    nobs      <- length(e)
    # pointer to the data matrix used.
    dataX     <- X
  }
  # END flag == TRUE

  # If data is not balanced:
  if (!flag){
    # IF using casewise selection:
    if (!pairwise){
      # If taking the balanced subset throws out all the data, quit.
      if (BnCS == 0 | BnTS == 0){
        stop("Either the number of CS observations per panel ",
              "or the number of TS observations per panel ",
              "used to compute the vcov matrix is zero. You must use",
              " pairwise selection.")
      }

      # Otherwise, estimate Sigma.hat using whole data and get
      # middle matrix (X' omega X).

      e <- balanced$resid        # residuals
      # Reshape the vector of residuals into the nCS x nTS matrix
      E <- matrix(e, BnCS, BnTS, byrow=TRUE)
      # Actually want it to be nTS x nCS
      E <- t(E)

      # Sigma.hat = (E'E)/BnTS
      Sigma.hat <- crossprod(E) / BnTS

      # If the average number of obs per panel used in the computation of the
      # pcses is less than half the avg number of obs in the original data,
      # suggest that the user use pairwise selection.
      if (avgN/2 > BnTS){
        warning("Caution! The number of CS observations per panel, ", BnTS,
                    ", used to compute the vcov matrix is less than half the",
                    "average number of obs per panel in the original data.",
                    "You should consider using pairwise selection.")
      }

      # X = rectangular data without first 3 cols (CS Id, TS Id, residuals)
      X           <- as.matrix(rect[ , 4:dim(rect)[2]])
      # If there are missings, make them 0 so R doesn't throw up.
      X[is.na(X)] <- 0
      # omega = Sigma.hat kronecker Identity_nTS
      omega       <- kronecker(Sigma.hat, diag(1, nTS))
      # middle = X' omega X
      middle      <- t(X) %*% omega %*% X
      # nobs = number of residuals in the original model
      nobs        <- length(resid(object))
    }
    # END Pairwise == FALSE

    # IF using pairwise selection:
    if (pairwise){
      # The next section of code follows gauss procedure from Franzeze (1996).

      # Get vector of 1/0 for valid row of obs or not.
      V <- rect[ , 4:dim(rect)[2]] # X data
      # vector = 1 if row is valid, 0 else.
      valid <- apply(!is.na.data.frame(V), 1, prod)
      nobs  <- sum(valid)          # number of obs = number of valid obs

      # Reshape valid and E to RnCS x RnTS
      e           <- rect$resid    # residuals
      e[is.na(e)] <- 0             # replace missing with 0
      E <- matrix(e, RnCS, RnTS, byrow=TRUE)
      # Actually, want RnTS x RnCS
      E <- t(E)

      V <- matrix(valid, RnCS, RnTS, byrow=TRUE)
      # Want RnTS x RnCS
      V <- t(V)

      # Numerator for Sigma.hat = E'E
      numer <- crossprod(E)
      # Denominator is number of observations per panel, V'V does this
      denom <- crossprod(V)
      # If not valid observation, replace with NA.
      denom[denom == 0] <- NA

      # If there is a set with no obs in common with other units, quit.
      check <- is.na(denom)
      if (sum(check) != 0){
        stop("Error! A CS-unit exists without any obs or without any obs in
              common with another CS-unit. You must remove that unit from the
              data passed to pcse().")
      }

      # Element by element division to get right T in denom.

      # Sigma.hat = (E'E) / (V'V) = (E'E) / T_j
      Sigma.hat <- numer/denom
      # X = rectangular data without first 3 cols (CS Id, TS Id, residuals)
      X           <- as.matrix(rect[ , 4:dim(rect)[2]])
      # If there are missings, make them 0 so R doesn't throw up.
      X[is.na(X)] <- 0
      # Omega = Sigma.hat kronecker Identity_T
      omega       <- kronecker(Sigma.hat, diag(1, nTS))
      # middle matrix = X' omega X
      middle      <- t(X) %*% omega %*% X
    }
    # End pairwise == TRUE
  }
  # End FLAG == FALSE

  # Finish up and print results.
  
  # estimate X'X with whole X.
  # estimate vcov with middle and X'X
  # calculate results.

  # X'X
  XX     <- t(X) %*% X
  # (X'X)^(-1)
  XXinv  <- solve(XX)
  # Variance covariance matrix = (X'X)^(-1) X' omega X (X'X)^(-1)
  vcov   <- XXinv %*% middle %*% XXinv
  # pcse = square root of diagonal elements of variance-covariance matrix
  pcse   <- sqrt(diag(vcov))
   
  # beta = output betas from 'lm' object
  b      <- summary(object)$coef[ , 1]
  # t-statistics = beta/pcse
  tstats <- b/pcse
  # degrees of freedom = number of observations - params in model
  # Note: cols of X = params in model, since X is model.matrix
  df     <- nobs - ncol(X)
  # p-value = 2*t-distribution with absolute val of t-stats for quantiles,
  #           degrees of freedom = df, no lower tail (want P[X > x])
  pval   <- 2*pt(abs(tstats), df, lower.tail=FALSE)

  # Bind results into a list.
  res <- list(vcov=vcov, pcse=pcse, b=b, tstats=tstats, df=df, pval=pval,
              pairwise=pairwise, nobs=nobs, nmiss=(nCS*nTS)-nobs, call=mc)

  # Attach class 'pcse' to results for summary method
  class(res) <- "pcse"
  # Return results.
  return(res)
}
# End pcse.
