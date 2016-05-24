autopls <- function (formula, data, testset = NULL, tselect = 'none', 
  prep = 'none', val = "LOO", scaling = TRUE, stingy = TRUE, 
  verbose = TRUE, backselect = 'auto', jt.thresh = 0.1,
  vip.thresh = 0.2, jump = NA, lower = NA, method = 'oscorespls')
{

  # Get the data matrices
  if (missing (data)) data <- environment (formula)
  mf <- match.call (expand.dots = FALSE)
  m <- match (c('formula', 'data'), names (mf), 0)
  mf <- mf [c(1, m)]
  mf[[1]] <- as.name ('model.frame')
  mf <- eval (mf, parent.frame ())
  mt <- attr (mf, 'terms')  
  Y <- model.response (mf, 'numeric')
  X <- model.matrix (mt, mf) [, -1]
  
  # Match arguments
  method <- match.arg (method, c('kernelpls', 'oscorespls'))
  tselect <- match.arg (tselect, c('none', 'active', 'passive'))
  val <- match.arg (val, c('LOO', 'CV'), several.ok = TRUE)
  prep <- match.arg (prep, c('none', 'bn', 'msc'))
  bsl <- match.arg (backselect, c('auto', 'no', 'A1', 'A2', 'A3', 
                    'A4', 'A5', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 
                    'C1'), several.ok = TRUE)
  backselect <- ifelse ('no' %in% bsl, FALSE, TRUE)
                           
  # Testset
  test <- ifelse (!is.null (testset), TRUE, FALSE) 
  
  if (test) # if there is a test set 
  {
    Xt <- X [testset,]
    X <- X [-testset,]    
    Yt <- Y [testset]
    Y <- Y [-testset]
  }
  
  # Check data suitability
  if (min (apply (X, 2, var)) == 0) stop ('Predictor(s) with no variance')
  if (min (apply (X, 1, var)) == 0) warning ('Observation(s) with no variance')

  # Determine acceptable number of latent vectors
  maxnlv <- ceiling (nrow (X) * 0.1)
  # With large matrices, reduce mxnlv:
  if (maxnlv > 5) maxnlv <- log (nrow (X), base = 2)

  # --- FUNCTIONS ---------------------------------------------------------- #

  # Determine number of latent vectors where RMSEval minimizes (nlv)
  getbest <- function (v)
  {

    # Truncate RMSE to a max. of 30 LV
    if (length (v) > 30) v <- v [1:30]
    # If RMSE are long enough: filter
    if (length (v) > 4)
    {
      criterion <- 0.03 * (max (v [1:4]) - min (v))
      # Low-pass filter (avoiding conflicts with package signal)
      lp <- stats::filter (v, c(1, 2, 1))
      # Slope
      slope <- lp [2:length (lp)] - lp [1:(length (lp) - 1)]      
      # if error decreases (almost) monotonically 
      if (max (slope, na.rm = TRUE) < criterion)
      {
        logic <- v > min (v) + criterion
        aopt <- match (FALSE, logic)
      }  
      # if error has a minimum
      else
      {
        # First minimum
        minflt <- match (FALSE, diff (lp) < 0)
        # Step back until 'punish factor' is met if there is any backw. increase
        act <- v > (v [minflt] + criterion)
        # Last value BELOW threshold before minflt
        if (sum (act [1:minflt]) > 0) 
          aopt <- minflt - match (TRUE, act [minflt:1]) + 2
        # If no backward increase take minfilt
        else aopt <- minflt
      }
    }
    # if rmse is too short for filtering
    else
    {
      criterion <- 0.03 * (max (v) - min (v))
      dif <- v [2:length (v)] + criterion - v [1:length (v) - 1]
      if (length (dif) == 1 & dif [1] < 0) aopt <- 2
      else aopt <- match (FALSE, dif < 0)
    }

    
    return (aopt)
  }

  # Validation to be run in function tryaround
  checkout <- function (selcode, cleanX, cleantestX)
  {      
    
    # Reduce the original data according to the selection
    Xchk <- cleanX [,selcode == 1]
    
    # Optional preprocessing
    if (prep != 'none') Xchk <- prepro (Xchk, method = prep)

    YXchk <- data.frame (Y = Y, X = I (Xchk))            
    
    # Run plsr 
    # ... with cross-validation
    if (tselect == 'none' | tselect == 'passive') 
      chk <- plsr (Y ~ X, data = YXchk, scale = scaling, 
        validation = 'LOO', method = method)              
    # ... without cross-validation only if tselect == 'active'
    else 
      chk <- plsr (Y ~ X, data = YXchk, scale = scaling, 
        method = method)      

    # Test set validation
    if (tselect != 'none')
    {  
      # Prepare test set
      Xchkt <- cleantestX [,selcode == 1]
      if (prep != 'none') Xchkt <- prepro (Xchkt, method = prep)
      cdv <- data.frame (Y = Yt, X = I (Xchkt))      
    
      # Derive nlv from test set results if tselect is 'active'      
      if (tselect == 'active')
      { 
        chk.rmse <- RMSEP (chk, estimate = 'test', newdata = cdv, 
          intercept = FALSE)$val [1,,] # test set errors for error estimate
        chk.best <- getbest (chk.rmse) # test set for nlv
      }
      # Derive nlv from CV if tselect is 'passive'
      else
      {
        chk.rmse <- RMSEP (chk, estimate = c('CV', 'test'), newdata = cdv, 
          intercept = FALSE)$val [1:2,,] # returns CV and test set errors    
        chk.best <- getbest (chk.rmse [1,]) # take CV for nlv
        chk.rmse <- chk.rmse [2,] # take test set for error estimate
      }
    }
    else # case of no test set 
    {
      chk.rmse <- RMSEP (chk, estimate = c('CV'),  
          intercept = FALSE)$val [1,,] # returns CV errors    
      chk.best <- getbest (chk.rmse)
    }
      
    # test statistics 
    # returns (1.) error at maxnlv, (2.) nlv where error minimizes and
    # (3.) error at that point
    if (length (chk.rmse) >= maxnlv)
      chk.val <- c(chk.rmse [maxnlv], chk.best, chk.rmse [chk.best])    
    
    # If the number of available LV's is smaller than the critical value
    # use the lowest absolute error instead of the error at maxnlv
    else 
      chk.val <- c(chk.rmse [chk.best], chk.best, chk.rmse [chk.best])           
    return (chk.val)
  }

  # Reduce autocorrelation
  ac <- function (vec, selcode, pow = 1)
  {
    # Function for finding local maxima
    locmax <- function (vec)
    {
      kern <- 5
      stp <- (kern - 1) / 2
      rca <- as.vector (vec)
      lthrca <- length (rca)

      # Treat positive and negative values separately
      rcapos <- rca
      rcaneg <- rca
      rcapos [rca <= 0] <- NA
      rcaneg [rca >= 0] <- NA
      rcaneg <- abs (rcaneg)

      # Function
      stepthrough <- function (vec, stp, lthrca)
      {
        loc2 <- vector ()
        for (i in (stp+1):(lthrca-stp))
        {
          getit <- vec [i] == max (vec [(i-stp):(i+stp)], na.rm = TRUE)
          loc2 <- c (loc2, getit)
        }
        # Deal with the ends
        loc1 <- vector ()
        for (j in 1:stp)
        {
          getit <- vec [j] == max (vec [1:(j+stp)], na.rm = TRUE)
          loc1 <- c (loc1, getit)
        }
        loc3 <- vector ()
        for (k in (lthrca-stp+1):lthrca)
        {
          getit <- vec [k] == max (vec [(k-stp):lthrca], na.rm = TRUE)
          loc3 <- c (loc3, getit)
        }
        loc <- c(loc1, loc2, loc3)
        loc <- which (loc == TRUE)
        return (loc)
      }
      locpos <- suppressWarnings (stepthrough (rcapos, stp, lthrca))
      locneg <- suppressWarnings (stepthrough (rcaneg, stp, lthrca))

      # Local maxima:
      loc <- sort (c(locpos, locneg))
      return (loc)
    }

    # Compute correlation matrix of predictors
    cmat <- cor (sX)

    # Mean nearest neighbor
    diag (cmat) <- -1
    nn <- apply (cmat, 1, max)
    mnn <- median (nn)  # Change this if another threshold is needed

    vec.thinout <- vec
    vec.thinout [selcode == 0] <- NA
    thin <- vector ()
    prog <- TRUE
    while (prog == TRUE)
    {
      # Get local maxima from regression coefficients
      loc <- locmax (vec.thinout)
      # Add this selection to overall selection
      thin <- c (thin, loc)
      # Correlation to local maxima larger than criterion?
      snn <- cmat [loc,] >= (mnn / pow)
      # Remove such predictors
      if (length (loc) > 1) vec.thinout [colSums (snn) > 0] <- NA
      else vec.thinout [snn == TRUE] <- NA
      # Remove previous local maxima
      vec.thinout [thin] <- NA
      # Check if there are at least two predictors left
      if (length (which (!is.na (vec.thinout))) < 3) prog <- FALSE
    }

    newvec <- rep (0, length (vec))
    newvec [thin] <- 1
    return (newvec)
  }

  # Dynamic significance filter
  dynp <- function (stp)
  {
    stp <- floor (stp)
    vec <- sort (jt)
    idx <- 1:length (vec)
    return (jt <= vec [idx == stp])  
  }
  
  tryaround <- function ()
  {
    
    # Functions used in tryaround ()
    A0 <- function ()
    {
      selmat <- matrix (1, nrow = 1, ncol = seln)
      selmat <- dynp (jump)        
      result <- list (selmat, 'A0')
      return (result)
    }
    
    A1 <- function (v) {v [jt >= jt.thresh] <- 0; v}
    A2 <- function (v) {v [vip < vip.thresh] <- 0; v}
    A3 <- function (v) {v [vip - jt < comb.thresh] <- 0; v}
    A4 <- function (v) {v <- dynp (seln * 0.9); v}
    A5 <- function (v) {v <- dynp (seln * 0.75); v}
    B1 <- function (v) {v <- ac (reg.coef, A1(v)); v}
    B2 <- function (v) {v <- ac (jt, A1(v)); v}
    B3 <- function (v) {v <- ac (vip, A1(v)); v}
    B4 <- function (v) {v <- ac (reg.coef, A2(v)); v}
    B5 <- function (v) {v <- ac (jt, A2(v)); v}
    B6 <- function (v) {v <- ac (vip, A2(v)); v}
    C1 <- function (v) {v <- ac (reg.coef, rep (1, seln), pow = 2); v}
        
    buildselmat <- function (compn)
    {
      selmat <- matrix (1, nrow = length (compn), ncol = seln)
      rownames (selmat) <- compn
      for (i in 1:length (compn))
        selmat [i,] <- do.call (compn [i], list (selmat [i,]))
      return (selmat)   
    }

    # Get factory fresh X and Xt
    if (counter == 1) 
    {
      sX.fresh <- X
      if (test) sXt.fresh <- Xt
      else sXt.fresh <- NULL
    }  
    else 
    {
      sX.fresh <- X [,prev.selection]
      if (test) sXt.fresh <- Xt [,prev.selection]
      else sXt.fresh <- NULL
    }
     
    seln <- ncol (sX.fresh)                                                            
            
    # Thresholds (change if appropriate)
    jt.thresh <- jt.thresh 
    vip.thresh <- vip.thresh # Note: VIP is scaled
    comb.thresh <- - 0.1

    if (!is.na (jump) & counter == 1 & seln > jump) jumpmode <- TRUE
    else jumpmode <- FALSE

    if (jumpmode) result <- A0()
    else
    { 
      if (!'auto' %in% bsl) selmat <- buildselmat (bsl)
      else selmat <- buildselmat (c('A1','A4'))
      compn <- rownames (selmat)

      # Remove solutions with < 2 remaining vars or without variable reduction
      sums <- rowSums (selmat)
      use <- (sums != ncol (selmat)) == (sums > 1)            
      selmat <- selmat [use,]      
      compn <- compn [use]
      
      if (sum (use) > 0)  
      {
        # Single solutions
        if (is.null (dim (selmat)))
        { 
          unq <- 1
          checkmat <- matrix (NA, nrow = 1, ncol = 3)
          checkmat <- checkout (selmat, sX.fresh, sXt.fresh)       
          # errors at maxnlv
          errors <- checkmat [1]    
          # nlv where errors in validation minimize
          comps <- checkmat [2]    
          # errors at minimum
          errors.min <- checkmat [3]
          result <- list (selmat, names (use) [use], errors.min)        
        }
        else 
        {
          code <- apply (selmat, 1, FUN = function(x) paste (x, collapse = ''))
          grps <- as.numeric (as.factor (code))        
          selmat <- selmat [order (grps),]
          grps <- sort (grps)
          unq <- match (1:max(grps), grps)
          
          # Matrix for the results of checkout ()
          checkmat <- matrix (NA, nrow = nrow (selmat), ncol = 3)
          rownames (checkmat) <- rownames (selmat)
          
          # Use function checkout only if selcode is met for the first time
          for (i in 1:nrow(checkmat))
          {
            if (i %in% unq) checkmat [i,] <- checkout (selmat [i,], 
              sX.fresh, sXt.fresh)
            else checkmat [i,] <- checkmat [i-1,]
          }
          
          # Rearrange checkmat and selmat
          checkmat <- checkmat [order (rownames (checkmat)),]
          selmat <- selmat [order (rownames (selmat)),]
          
          # errors at maxnlv
          errors <- checkmat [,1]    
          # nlv where errors in validation minimize
          comps <- checkmat [,2]    
          # errors at minimum
          errors.min <- checkmat [,3]
          
          if (stingy)
          {
            # Scale errors
            scalederrors <- errors / sd (Y)
            # Select predictor sets leading to a scaled error that does not exceed
            # the minimum scaled error + x
            lowerrors <- scalederrors <= min (scalederrors) + 0.025
            # Remove models with higher error
            comps [!lowerrors] <- 9999
            # Get the models with the lowest nlv where error minimizes
            these <- comps == min (comps)
            scalederrors [!these] <- 9999
            # Among these models, get the first one with the lowest error
            thisone <- which.min (scalederrors)
          }
          else # several solutions and not stingy
          {
            # Get the models with the lowest error
            these <- errors.min == min (errors.min)          
            # Among these models, get the one with the lowest nlv
            comps [!these] <- 9999
            thisone <- which.min (comps)
          }        
          result <- list (selmat [thisone,], compn [thisone], 
            errors.min [thisone])
        }            
      }
      # Case of no useful solution
      else result <- NULL      
    }
    return (result)
  }

  # VIP.R: Implementation of VIP (variable importance in projection)(*) for the
  # `pls' package.
  # Copyright 2006,2007 Bjorn-Helge Mevik
  # This program is free software; you can redistribute it and/or modify
  # it under the terms of the GNU General Public License version 2 as
  # published by the Free Software Foundation.

  # (*) As described in Chong, Il-Gyo & Jun, Chi-Hyuck, 2005, Performance of
  # some variable selection methods when multicollinearity is present,
  # Chemometrics and Intelligent Laboratory Systems 78, 103--112.

  VIP <- function (object)
  {
      SS <- c(object$Yloadings) ^ 2 * colSums (object$scores^2)
      Wnorm2 <- colSums (object$loading.weights ^ 2)
      SSW <- sweep (object$loading.weights ^ 2, 2, SS / Wnorm2, "*")
      sqrt (nrow (SSW) * apply (SSW, 1, cumsum) / cumsum (SS))
  }  

  # -- END FUNCTIONS ------------------------------------------------------- #

  # At this point, a loop begins that stores model results after predictor
  # selection and turns back here.

  counter <- 1
  loop <- TRUE
  onemoreloop <- FALSE
  
  while (loop == TRUE | onemoreloop == TRUE)
  {

    # --- Variables from previous run and update of X ---------------------- #

    if (counter > 1)
    {
      # Previous selection of predictors
      prev.selection <-
        as.logical (unlist (tmp [[6 + counter]][length (model) + 5]))
      # Previous RMSE
      prev.rmse <- tmp [[6 + counter]] [[length (model) + 1]]
      # Previous R2
      prev.r2 <- tmp [[6 + counter]] [[length (model) + 2]]
      # Previous nlv
      prev.nlv <- tmp [[6 + counter]] [[length (model) + 3]]
      # Previous croperror
      prev.croperror <- tmp [[6 + counter]] [[length (model) + 4]]
      # Construct new X
      sX <- X [,prev.selection] #selection
      # Construct new X for the test set
      if (test) sXt <- Xt [,prev.selection]    
    }
    else # if in the first run rename X and Xt nevertheless
    {
      sX <- X 
      if (test) sXt <- Xt
    }
    
    # --- Data preparation ------------------------------------------------  #

    # X preprocessing
    if (prep != 'none') sX <- prepro (sX, method = prep) 
    # Build the dataframe
    sYX <- data.frame (Y = Y, X = I (sX))
    
    # --- Run PLSR --------------------------------------------------------- #
    
    model <- plsr (Y ~ X, data = sYX, jackknife = TRUE, 
      scale = scaling, validation = val, method = method)

    # Model parameters
    if (!test) # No test set
    {
      # Call the RMSEP and R2 methods for mvr objects
      rmse <- RMSEP (model, estimate = c('train', 'CV'), 
        intercept = FALSE) $val [1:2,,]
      r2 <- R2 (model, estimate = c('train', 'CV'), 
        intercept = FALSE) $val [1:2,,]
      nlv <- getbest (rmse [2,])
      rmse.cal <- rmse [1,nlv]
      rmse.val <- rmse [2,nlv]
      r2.cal <- r2 [1,nlv]
      r2.val <- r2 [2,nlv]
    }
    else # Test set present
    {                                                          
      # Build test set data and call methods for class mvr      
      if (prep != 'none') sXt <- prepro (sXt, method = prep)
      valdat <- data.frame (Y = Yt, X = I (sXt))
      rmse <- RMSEP (model, estimate = c('train', 'CV', 'test'), 
        newdata = valdat, intercept = FALSE) $val [1:3,,]
      r2 <- R2 (model, estimate = c('train', 'CV', 'test'), 
        newdata = valdat, intercept = FALSE) $val [1:3,,]
      
      # selected nlv from test or CV
      if (tselect == 'active') nlv <- getbest (rmse [3,]) # from test set
      else nlv <- getbest (rmse [2,]) # from cross-validation 
      
      rmse.cal <- rmse [1,nlv]
      rmse.val <- rmse [2,nlv]
      rmse.test <- rmse [3,nlv]
      r2.cal <- r2 [1,nlv]
      r2.val <- r2 [2,nlv]
      r2.test <- r2 [3,nlv]
    }
    
    reg.coef <- as.vector (coef (model, ncomp = nlv))
    jt <- as.vector (suppressWarnings (jack.test (model, ncomp = nlv) $pvalues))
    vip <- VIP (model) [nlv,]
    # Rescale VIP
    vip <- (vip - min (vip)) / max (vip - min (vip))
    
    # ---- Reporting ------ #
    
    if (verbose)
    {
      if (counter == 1) insp <- nchar (ncol (X)) + 2
      spaces <- 
        function (sp, item) paste (rep (' ', sp - nchar (item)), collapse = '')
      space1 <- spaces (4, counter)
      space2 <- spaces (insp, ncol (sX))
      space3 <- spaces (4, nlv)
      space4 <- spaces (7, round (r2.val, 3))
      space6 <- spaces (7, round (rmse.val, 3))
      if (test) 
      {
        space5 <- spaces (7, round (r2.test, 3))
        space7 <- spaces (7, round (rmse.test, 3))
      }
      
      cat (counter, space1, 
           'Pred: ', ncol (sX), space2, 
           'LV: ', nlv, space3,
           'R2v: ', round (r2.val, 3), space4, 
           sep = '')
      if (test) cat ('R2t: ', round (r2.test, 3), space5, sep = '')
      cat ('RMSEv: ', round (rmse.val, 3), space6, sep = '')
      if (test) cat ('RMSEt: ', round (rmse.test, 3), space7, sep = '')
      if (counter > 1) cat (paste ('Criterion:',  res [[2]]))
      cat ('\n')
      flush.console ()
    }
     
    # ---------- Selection of predictors for subsequent runs --------------- #

    selcode <- NA

    if (backselect & !onemoreloop)
    {
      res <- tryaround ()
      if (!is.null (res)) selcode <- unlist (res [[1]])    
      else backselect <- FALSE
    }
    
    if (counter > 1)
    {
      new.selection <- prev.selection      
      new.selection [prev.selection == TRUE] <- selcode
      selcode <- new.selection
    }
    
    # --- Result object ---------------------------------------------------- #

    # --- Construct tables of RMSE, R2 and a vector of nlv
    # --- These values are used in the final selection of an iteration

    # Construct rmse.crop for use with stingy
    if (nlv > maxnlv) # case of large nlv
    {
      # If tselect is not 'none' use test set errors for the final selection
      #
      if (tselect != 'none') 
      {
        rmse.crop <- rmse [3, maxnlv] # in test set
        r2.crop <- r2 [3, maxnlv] # in test set
      }
      else # otherwise use CV
      {
        rmse.crop <- rmse [2, maxnlv] # in CV
        r2.crop <- r2 [2, maxnlv] # in CV
      }
    }    
    # In case of small nlv
    else 
    {
      if (tselect != 'none') 
      {
        rmse.crop <- rmse [3, nlv] # in test set
        r2.crop <- r2 [3, nlv] # in test set
      }
      else
      {
        rmse.crop <- rmse [2, nlv] # in CV
        r2.crop <- r2 [2, nlv] # in CV
      }
    }
    
    if (test)
    {
      metarmse <- c(rmse.cal, rmse.val, rmse.crop, rmse.test)
      names (metarmse) <- c('RMSEcal', 'RMSEval', 'RMSEcrop', 'RMSEtest')
      metar2 <- c(r2.cal, r2.val, r2.crop, r2.test)
      names (metar2) <- c('R2cal', 'R2val', 'R2crop', 'R2test')
    }
    else
    {
      metarmse <- c(rmse.cal, rmse.val, rmse.crop)
      names (metarmse) <- c('RMSEcal', 'RMSEval', 'RMSEcrop')
      metar2 <- c(r2.cal, r2.val, r2.crop)
      names (metar2) <- c('R2cal', 'R2val', 'R2crop')
    }
     
    metanlv <- nlv

    if (counter > 1)
    {
      # append RMSE to the existing table of RMSE's
      iter.1 <- prev.rmse
      metarmse <- cbind (iter.1, metarmse)
      colnames (metarmse) [counter] <- paste ('iter.', counter, sep = '')

      # same for R2
      iter.1 <- prev.r2
      metar2 <- cbind (iter.1, metar2)
      colnames (metar2) [counter] <- paste ('iter.', counter, sep = '')

      # same for nlv
      iter.1 <- prev.nlv
      metanlv <- c (iter.1, nlv)
      names (metanlv) [counter] <- paste ('iter.', counter, sep = '')
    }

    # --- First time here and in backselection mode? Prepare result object

    if (counter == 1)
    {
      # Collect information about the models (future 'metapls' object)
      tmp <- list ('preprocessing' = prep,       # [[1]]
                   'scaling'       = scaling,    # [[2]]
                   'RMSE'          = NA,         # [[3]]
                   'R2'            = NA,         # [[4]]
                   'nlv'           = NA,         # [[5]]
                   'best.model'    = NA,         # [[6]]
                   'predictors'    = NA)         # [[7]]

      predictors <- rep (TRUE, ncol (sX)) # All predictors have been used
      names (metanlv) [1] <- 'iter.1'
    }
    else predictors <- prev.selection

    # List element containing the iteration results

    index <- 7 + counter
    lth <- length (model)

    tmp [[index]] <- model
    names (tmp)[index] <- paste ('iter.', counter, sep = '')

    # Element [[1]]: RMSE
    tmp [[index]][[lth + 1]] <- metarmse
    names (tmp [[index]])[lth + 1] <- 'RMSE'

    # Element [[2]]: R2
    tmp [[index]][[lth + 2]] <- metar2
    names (tmp [[index]])[lth + 2] <- 'R2'

    # Element [[3]]: Number of latent vectors (nlv)
    tmp [[index]][[lth + 3]] <- metanlv
    names (tmp [[index]])[lth + 3] <- 'nlv'

    # Element [[4]]: Predictor selection for current run
    tmp [[index]][[lth + 4]] <- predictors
    names (tmp [[index]])[lth + 4] <- 'predictors'
    
    # Element [[5]]: Predictor selection for next run
    tmp [[index]][[lth + 5]] <- selcode
    names (tmp [[index]])[lth + 5] <- 'selection'
    
    # Continue?
    if (onemoreloop) 
    {
      loop <- FALSE
      onemoreloop <- FALSE
    }  
    else
    {
      if (counter > 1)
      {
        lth <- ncol (metar2)      
        loop <- FALSE
        
        # Continue if the r2 in the last model is not worse than "limit"
        # In case of tselect != 'none', this refers to the test set validation
        
        if (tselect != 'none') 
        {
          if (!is.na (lower)) limit <- lower 
          else limit <- max (metar2 [4,]) - 0.05 # in test set
          if (limit > 0 & metar2 [4,lth] > limit) loop <- TRUE
        }  
        else 
        {
          if (!is.na (lower)) limit <- lower 
          else limit <- max (metar2 [2,]) - 0.05 # in CV
          if (limit > 0 & metar2 [2,lth] > limit) loop <- TRUE
        }
        # Continue also if the nlv decreases    
        if (metanlv [lth] < metanlv [lth - 1]) loop <- TRUE
      }

      # exit if backselect = 'no'
      if (!backselect) loop <- FALSE
      
      # action if less than 10 predictors have been selected
      else if (sum (selcode) < 10) loop <- FALSE
          
      # Now, if loop == F and one of the tested models performed 
      # better than previous models one more loop should take place but without 
      # further backwards selection.
      
      if (backselect & loop == FALSE)
      {
        if (tselect != 'none') 
        {
          if (counter > 1) limit <- max (metarmse [4,]) # in test set
          else  limit <- max (metarmse [4])
          if (limit > 0 & res [[3]] < limit) onemoreloop <- TRUE
        }  
        else 
        {
          if (counter > 1) limit <- max (metarmse [2,]) # in CV
          else limit <- max (metarmse [2]) # in CV        
          if (limit > 0 & res [[3]] < limit) onemoreloop <- TRUE
        }
      }
    }
    
    # Final action in case of no return
    if (loop == FALSE & onemoreloop == FALSE)
    {      
      
      # Parameters relating to the best model
      # Case of no selection but backselect == TRUE or numeric
      if (counter == 1)
      {
        best <- 1
        tmp [[7]] <- predictors                           # selection
        tmp [[6]] <- best                                 # best model
        tmp [[5]] <- tmp [[8]][length (model) + 3]        # nlv
        tmp [[4]] <- metar2                               # R2
        tmp [[3]] <- metarmse                             # RMSE
      }
      # Case of prev. selection and backselect == TRUE or numeric
      else
      {
        # ---- Final model selection out of all iterations:

        if (tselect != 'none')
        {
          meta.error <- unlist (metarmse [4,]) # based on test set at nlv
          meta.crop <- unlist (metarmse [3,]) # based on test set at maxnlv
        }
        else
        {
          meta.error <- unlist (metarmse [2,]) # based on CV at nlv
          meta.crop <- unlist (metarmse [3,]) # based on CV at maxnlv
        }


        # this is the first model with the global, minimum RMSEval:
        errorsel <- which (meta.error == min (meta.error))[1]

        if (stingy)
        {
          # nlv in this model:
          errorsel.nlv <- tmp [[7 + errorsel]][[length (model) + 3]][errorsel]
          # If the model with the lowest RMSEval has a nlv <= maxnlv take that:
          if (errorsel.nlv <= maxnlv) best <- errorsel
          # otherwise take the model with the lowest RMSE at maxnlv:
          else best <- which (meta.crop == min (meta.crop))[1]
        }
        else best <- errorsel
        names (best) <- NULL
      }
  
      # Model parameters
      opt.lv <- metanlv [best] # sequence of lv
      names (opt.lv) <- NULL
      preds <- unlist (tmp [[7 + best]][length (model) + 4]) # last predictors
      names (preds) <- NULL      
      
      tmp [[7]] <- preds # last predictors
      tmp [[6]] <- best  # selected model
      tmp [[5]] <- opt.lv # nlv

      if (is.null (dim (metarmse) [2]))
      {
        opt.rmse <- metarmse [2]
        opt.r2 <- metar2 [2]
        tmp [[4]] <- opt.r2
        tmp [[3]] <- opt.rmse
      }
      else
      {
        opt.rmse <- metarmse [2, best]
        opt.r2 <- metar2 [2, best]
        tmp [[4]] <- metar2 [,best]  
        tmp [[3]] <- metarmse  [,best] 
      }

      # The result object consists of three parts: 
      # (1) The first part represents the "best" model object and has the 
      #     structure of an mvr object in package pls.
      # (2) The second is a list object called metapls and contains desriptive
      #     statistics based on the selected number of latent vectors,
      #     predictors used, summarizing statistics from the other 
      #     iterations and the entire, unreduced set of predictors
      # (3) The third consists of other iterations and can be used to produce 
      #     a new autopls object based on another iteration.
      
      # Get "best" model and put it in part 1
      part1 <- tmp [[7 + best]][!as.logical(match(names(tmp [[7 + best]]),
          c("RMSE", "R2", "nlv", "selection"),0))]
      
      # Construct part 2 ($metapls)  
      part2 <- list (current.iter  = best,
                     autopls.iter  = best,
                     current.lv    = opt.lv,       
                     autopls.lv    = opt.lv,
                     lv.history    = tmp [[length (tmp)]] $nlv,
                     rmse.history  = tmp [[length (tmp)]] $RMSE,
                     r2.history    = tmp [[length (tmp)]] $R2,
                     X             = X,
                     Y             = Y,
                     X.testset     = NULL, 
                     Y.testset     = NULL,
                     preprocessing = prep,
                     scaling       = scaling,
                     val           = val,
                     call          = sys.call ())
      
      # Add test set if applicable
      if (test) 
      {
        part2$X.testset <- Xt
        part2$Y.testset <- Yt
      }  
      
      # Get data for part 3 ($iterations)
      part3 <- tmp [8:length(tmp)]
      
      # In order to save memory, remove original data from all iterations
      # This can be reconstructed using $X and $Y ($metapls) and $predictors
      for (i in 1:length(part3))
      { 
        part3 [[i]] <- part3 [[i]] [!as.logical (match (names 
          (tmp [[7 + best]]), c("RMSE", "R2", "nlv", "selection"),0))]
        part3 [[i]] $model$Y <- NULL
        part3 [[i]] $model$sX <- NULL
      }
      
      # Put parts together
      result <- c(part1)
      result [[length (result) + 1]] <- part2
      names (result) [length (result)] <- 'metapls'
      result [[length (result) + 1]] <- part3
      names (result) [length (result)] <- 'iterations'
    }
    
    # Update counter
    counter <- counter + 1
  }
  
  class (result) <- 'autopls'

  # Reporting
  if (verbose) print (result)  
  invisible (result)
}
