############################################################################
######################### bcbcsf main functions ############################
############################################################################
bcbcsf_fitpred <- function (
  ## arguments specifying info of data sets
  X_tr, y_tr, nos_fsel = ncol (X_tr), 
  X_ts = NULL,  standardize = FALSE, rankf = FALSE,
  ## arguments for prediction
  burn = NULL, thin = 1, offset_sdxj = 0.5,
  ## arguments for Markov chain sampling
  no_rmc = 1000, no_imc = 5, no_mhwmux = 10,
  fit_bcbcsf_filepre = ".fitbcbcsf_", 
  ## arguments specifying priors for parameters and hyerparameters
  w0_mu = 0.05, alpha0_mu = 0.5, alpha1_mu = 3,
  w0_x  = 1.00, alpha0_x  = 0.5, alpha1_x  = 10,
  w0_nu = 0.05, alpha0_nu = 0.5, prior_psi = NULL,
  ## arguments for metropolis sampling for wmu, wx
  stepadj_mhwmux = 1, diag_mhwmux = FALSE,
  ## arguments for computing adjustment factor
  bcor = 1, cut_qf = exp (-10), cut_dpoi = exp (-10), nos_sim = 1000,
  ## whether look at progress
  monitor = TRUE)
{
  if (!is.matrix (X_tr)) stop ("X_tr must be a matrix")

  ## read information about data
  n <- nrow (X_tr) ## numbers of obs
  p <- ncol (X_tr)

  ## find number of observations in each group, posterior freq of y
  nos_g <- as.vector (table (y_tr))
  G <- length (nos_g)
  if (any(nos_g < 2)) stop ("less than 2 cases in some group in your data")
  
  ## set prior for proportion of cases
  if (is.null (prior_psi)) prior_psi <- rep (1, G)
  if (length (prior_psi) != G) stop ("length of prior_psi is wrong")

  ## prior frequency of class lables
  post_y <- nos_g + prior_psi
  freqy <- post_y / sum (post_y)

  ################### Preprocessing Features ##############################
  ## 1) feature selection on original data (the same as standardized data)
  info_sel <- list (vars = 1:p)
  if (any (c(rankf, nos_fsel < p)))
    info_sel <- rank_F (X_tr, y_tr)

  ## 2) ordering and selecting features
  p_fselmax <- max (nos_fsel)
  fselmax <- info_sel $ vars [1 : p_fselmax]
  X_tr <- X_tr[, fselmax, drop = FALSE]

  ## 3) standardizing retained features
  mle_ori_fselmax <- trpr_mle (
    X_tr = X_tr, y_tr = y_tr, rankf = FALSE)$list_fit_mle[[1]]

  nuj_ori_fselmax <- rep (0, p_fselmax)
  wxj_ori_fselmax <- rep (1, p_fselmax)
  if (standardize)
  {
    nuj_ori_fselmax <- mle_ori_fselmax $ nuj
    wxj_ori_fselmax <- mle_ori_fselmax $ wxj
    X_tr <- sweep (X_tr, 2, nuj_ori_fselmax, "-")
    X_tr <- sweep (X_tr, 2, sqrt (wxj_ori_fselmax), "/")
  }

  ## 4) Gathering sufficient statistic on standardized data
  tgsum_X_fselmax <- t (rowsum (X_tr,y_tr)) ## grouped sum of features
  sum_X2_fselmax <- colSums (X_tr^2)
  #########################################################################

  ## creating storage of results
  nnfsel <- length (nos_fsel)
  fitfiles <- rep ("", nnfsel)

  array_probs_pred <- NULL
  if (!is.null (X_ts))
  {
    n_ts <- nrow (X_ts)
    array_probs_pred <- array (0, dim = c(n_ts, G, nnfsel) )
  }

  if (monitor)  
  cat ("   Be Patient ... BCBCSF is fitting...  \n")
  
  finished <- 0
  total <- sum (nos_fsel) * no_rmc
  if (monitor)  pb <- txtProgressBar(min = 0,  max = total, style = 3)
  #########################################################################
  ## starting training and prediction for each number of retained features
  for (i  in seq (1, nnfsel) )
  {
      ## information on selected k features
      k <- nos_fsel [i]
    
      if (!is.null (fit_bcbcsf_filepre))
      fitfiles[i] <-  
          paste (fit_bcbcsf_filepre,
            "alpha1_mu_",alpha1_mu, "_n_",n, "_p_", p, 
            "_nfsel_", k, "_biascor_", bcor, ".RData", sep = "")

      ## Start Marlov chain super-transition
      if (k >= 1) 
      {
        ## information on feature selection and standardization
        fsel <- info_sel $ vars [1:k]
        nuj_std <- nuj_ori_fselmax [1:k]
        wxj_std <- wxj_ori_fselmax [1:k]
        cut_F <- info_sel $ fstats [k]
        nos_omit <- p - k

        ## sufficient statstic on selected features
        tgsum_X <- tgsum_X_fselmax [1:k,,drop = FALSE]
        sum_X2 <- sum_X2_fselmax [1:k]

        ## compute qf and partial lambda for bias correction
        if (bcor == 1 & k < p)
        {
          qflmd <- gen_qflmd (y_tr, cut_F, alpha1_mu, alpha1_x, 
                              cut_qf, nos_sim)
        }

        ## static variables in Gibbs sampling
        alpha_wmuj <- (alpha1_mu + G) / 2
        alpha_wxj <- (alpha1_x + n) / 2
        alpha_wnu <- (alpha0_nu + k) / 2
        alpha_wmu <- alpha1_mu * k / 2  - alpha0_mu / 2
        alpha_wx <- alpha1_x * k / 2  - alpha0_x / 2
        lambda0_wmu <- alpha0_mu * w0_mu / 2 ##lambda for wmu from prior
        lambda0_wx <- alpha0_x * w0_x / 2 ##lambda for wx from prior

        ## Metropolis Sampling stepsize
        stepsizes_mhwmux <- stepadj_mhwmux /
            sqrt ( 10 * c(max(alpha0_mu,alpha_wmu), max(alpha0_x, alpha_wx)) )

        ## initialize Markov chain from MLE
        muj    <- mle_ori_fselmax$muj[1:k,,drop = FALSE] - nuj_std
        wxj    <- mle_ori_fselmax$wxj[1:k] / wxj_std
        wmuj   <- mle_ori_fselmax$wmuj[1:k] * 0.01
        nuj    <- rowMeans (muj)
        wx     <- 1/mean (1/wxj)
        logwx  <- log (wx)
        wmu    <- w0_mu
        logwmu <- log (wmu)
        wnu    <- 1

        ## Markov chain storage
        MUJ   <- array (0, dim = c(k, G, no_rmc))
        WXJ   <- array (0, dim = c(k, no_rmc))
        WMUJ  <- array (0, dim = c(k, no_rmc))
        NUJ   <- array (0, dim = c(k, no_rmc))
        WX    <- array (0, dim = no_rmc) ## a vector
        WMU   <- array (0, dim = no_rmc) ## a vector
        WNU   <- array (0, dim = no_rmc) ## a vector

        ## start Gibbs sampling
        j_save <- 1
        i_save <- no_imc * j_save
        for (i_mc in 1 : (no_rmc * no_imc))
        {
          ## update muj
          vars_muj <- 1 / (1/wmuj + outer(1/wxj,nos_g) )
          means_muj <- (nuj / wmuj + tgsum_X / wxj) * vars_muj
          muj <- means_muj + replicate (G, rnorm (k)) * sqrt (vars_muj)

          ## update wxj
          lambda_wxj <-  (alpha1_x * wx + sum_X2 -
                          2 * rowSums (tgsum_X * muj) + 
                          rowSums (sweep (muj^2, 2, nos_g, "*")) )/2 

          wxj <- rinvgam (k, alpha_wxj, lambda_wxj)
          ## update wmuj
          lambda_wmuj <- (alpha1_mu * wmu + rowSums ((muj - nuj)^2) ) / 2
          wmuj <- rinvgam (k, alpha_wmuj, lambda_wmuj)

          ## update nuj
          sum_muj <- rowSums (muj)
          vars_nuj <- 1 / (1 / wnu + G / wmuj)
          means_nuj <- sum_muj / wmuj * vars_nuj
          nuj <- means_nuj + rnorm (k) * sqrt (vars_nuj)

          ## update wnu
          lambda_wnu <- (alpha0_nu * w0_nu + sum (nuj^2) ) / 2
          wnu <- rinvgam (1, alpha_wnu, lambda_wnu)

          ## update wx and wu together with M-H methods
          ## log posterior of log (wmu, wx)
          lambda_wmu <- alpha1_mu * sum(1/wmuj) / 2
          lambda_wx <- alpha1_x * sum(1/wxj) / 2
          logpost_logwmux <-  function (lw)
          {
            w <- exp (lw)
            b4cor <-
              alpha_wmu * lw[1] - lambda_wmu * w[1] -
              lambda0_wmu / w[1] + alpha_wx * lw[2] -
              lambda_wx * w[2] - lambda0_wx / w[2]

            if (bcor == 1 & nos_omit > 0)
              b4cor + nos_omit *
              log (comp_adjfactor (w[1], w[2], qflmd, cut_dpoi))
            else b4cor
          }

          log_wmu_wx <- met_gauss (
            iters = no_mhwmux, log_f = logpost_logwmux,
            ini_value = c(logwmu, logwx), stepsize = stepsizes_mhwmux,
            diag_mh = diag_mhwmux )
          logwmu <- log_wmu_wx [1]
          logwx <- log_wmu_wx [2]
          wmu <- exp (logwmu)
          wx <- exp (logwx)

          ## write states into Marlov chain arrays
          if (i_mc == i_save)
          {
            MUJ [,,j_save] <- muj
            WXJ [,j_save] <- wxj
            NUJ [,j_save] <- nuj
            WMUJ [,j_save] <- wmuj
            WX [j_save] <- wx
            WMU [j_save] <- wmu
            WNU [j_save] <- wnu
            j_save <- j_save + 1
            i_save <- j_save * no_imc
            finished <- finished + nos_fsel [i]
            if (monitor)   
            {
                setTxtProgressBar(pb, finished)
            }
          }
        }

        fit_bcbcsf <- list (
          fsel =  fsel, nuj_std = nuj_std, wxj_std = wxj_std,
          MUJ = MUJ, WXJ = WXJ, NUJ = NUJ, WMUJ = WMUJ, WX = WX,
          WMU = WMU, WNU = WNU, freqy = freqy,
          no_imc =  no_imc, no_rmc = no_rmc, bias_corrected = bcor )
      }
      else  fit_bcbcsf <- list (fsel = NULL, freqy = freqy)
      
      if (fitfiles[i] != "") save (fit_bcbcsf, file = fitfiles[i])

      ###################### making prediction ###############################
      if (!is.null (X_ts))
      {
          array_probs_pred[,,i] <- 
          mcmc_pred (X_ts = X_ts, fit_bcbcsf = fit_bcbcsf,
                     burn = burn, thin = thin, offset_sdxj = offset_sdxj)
      }
  }
  
  if (monitor) close (pb)

  if (!is.null (array_probs_pred))
  {
      dims <- dim (array_probs_pred)
      dimnames (array_probs_pred)  <- list(paste("Case", 1:dims[1], sep=""),
                                           paste("Class", 1:dims[2], sep=""),
                                           paste("fsel", 1:dims[3], sep=""))
  }
  
  ## returning results
  list (fit_bcbcsf = fit_bcbcsf,
        fitfiles = fitfiles,
        array_probs_pred = array_probs_pred,
        nos_fsel = nos_fsel )
}


############################################################################
######################### functions for prediction #########################
############################################################################

bcbcsf_pred <- function ( 
      X_ts, out_fit, burn = NULL, thin = 1, offset_sdxj = 0.5)
{
    if (is.vector (X_ts)) X_ts <- matrix (X_ts,1,)
    
    fitfiles <- out_fit$fitfiles
    nos_fsel <- out_fit$nos_fsel
    
    array_probs_pred <- NULL
    
    for (i in 1:length (nos_fsel))
    {
      fit_bcbcsf <- reload_fit_bcbcsf (fitfiles[i])
      probs_pred <- mcmc_pred (X_ts = X_ts, fit_bcbcsf = fit_bcbcsf,
                    burn = burn, thin = thin, offset_sdxj = offset_sdxj)
      array_probs_pred <- abind (array_probs_pred, probs_pred, along = 3)
    }
    
    dims <- dim (array_probs_pred)
    dimnames (array_probs_pred)  <- list(paste("Case", 1:dims[1], sep=""),
                                         paste("Class", 1:dims[2], sep=""),
                                         paste("fsel", 1:dims[3], sep=""))
    
    
    list (fitfiles = fitfiles,
          array_probs_pred = array_probs_pred, 
          nos_fsel = nos_fsel)
}

mcmc_pred <- function (
    X_ts, fit_bcbcsf = NULL, fit_bcbcsf_file = NULL,
    burn = NULL,  thin = 1, offset_sdxj = 0.5)
{
  n <- nrow (X_ts)

  if (is.null (fit_bcbcsf))
  {
    fit_bcbcsf <- reload_fit_bcbcsf (fit_bcbcsf_file)
  }

  if (is.null (fit_bcbcsf$fsel)) ## no features used
  {
    t (replicate (n, fit_bcbcsf$freqy) )
  }
  else
  {
    if (is.null (burn)) burn <- floor (fit_bcbcsf$no_rmc * 0.2)

    mu_dim <- dim (fit_bcbcsf$MUJ)
    k <- mu_dim [1]
    G <- mu_dim [2]
    no_rmc <- mu_dim [3]

    ## standardizing and selecting features
    X_ts <- X_ts[, fit_bcbcsf$fsel, drop = FALSE]
    X_ts <- sweep (X_ts,2, fit_bcbcsf$nuj_std, "-")
    X_ts <- sweep (X_ts,2, sqrt(fit_bcbcsf$wxj_std), "/")

    ## prepare indice of samples used to predict
    ix_pred <- burn + thin * seq(0, (no_rmc - burn) %/% thin )

    SDXJ <- sqrt(fit_bcbcsf$WXJ[,ix_pred, drop = FALSE])
    if (offset_sdxj > 1E-5)
    {
        offset <- quantile (SDXJ, offset_sdxj)
    }
    else offset <- 0
    
    SDXJ <- SDXJ + offset

    .C ( "pred_ht", n, k, G, length(ix_pred), X_ts,
          fit_bcbcsf$MUJ[,,ix_pred], SDXJ, log(fit_bcbcsf$freqy),
          probs_pred = matrix(0,n,G), PACKAGE = "BCBCSF"
      ) $ probs_pred
  }
}

mlepred <- function (X_ts, fit_mle)
{
  n <- nrow (X_ts)

  if (is.null (fit_mle$fsel) )
  {
    t(replicate (n, fit_mle$freqy))
  }
  else
  {
    k <- nrow (fit_mle$muj)
    G <- ncol (fit_mle$muj)

    .C ("pred_ht",
         n, k, G, as.integer(1), X_ts[, fit_mle$fsel], fit_mle$muj,
         sqrt (fit_mle$wxj),log (fit_mle$freqy),
         probs_pred = matrix (0, n, G), PACKAGE = "BCBCSF"
       ) $ probs_pred
  }
}

############################################################################
######################### BCBCSF result Analyzing ##########################
############################################################################


bcbcsf_sumfit <- function (
  fit_bcbcsf = NULL, fit_bcbcsf_afile = NULL, burn = NULL, thin = 1)
{

  if (is.null(fit_bcbcsf)) fit_bcbcsf <- reload_fit_bcbcsf (fit_bcbcsf_afile)

  if (is.null (burn)) burn <- floor (fit_bcbcsf$no_rmc * 0.2)

  nuj <- apply (fit_bcbcsf $ NUJ[, - (1:burn), drop = FALSE], 1, median )
  wmuj <- apply (fit_bcbcsf $ WMUJ[, - (1:burn), drop = FALSE], 1, median )
  wx <-  median (fit_bcbcsf $ WX[- (1:burn)] )
  wmu <-  median (fit_bcbcsf $ WMU[- (1:burn)] )
  wnu <-  median (fit_bcbcsf $ WNU[- (1:burn)] )

  muj <- apply (fit_bcbcsf $ MUJ[,, - (1:burn), drop = FALSE], c(1,2), median)
  wxj <- apply (fit_bcbcsf $ WXJ[,  - (1:burn), drop = FALSE], 1, median)
  cmuj <- muj - apply (muj,1, mean)
  scmuj <- cmuj/sqrt(wxj)
  signalj <- apply (scmuj, 1, sd)
  
  freqy <- fit_bcbcsf$freqy
  fsel <- fit_bcbcsf$fsel

  list (
    nuj_std = fit_bcbcsf$nuj_std, 
    wxj_std = fit_bcbcsf$wxj_std,
    nuj = nuj, 
    wx = wx, 
    wmu = wmu, 
    wnu = wnu,
    wmuj = wmuj,  
    cmuj = cmuj, 
    muj = muj, 
    wxj = wxj,  
    scmuj = scmuj,
    signalj = signalj,
    freqy = freqy, 
    fsel = fsel 
  )
}

reload_fit_bcbcsf <- function (fit_bcbcsf_afile)
{
    local ({
            fit_bcbcsf <- get(load (fit_bcbcsf_afile))
            return (fit_bcbcsf)
          })
}

bcbcsf_plotsumfit <- function (sum_fit)
{
   G <- ncol (sum_fit$scmuj)
   par (mfrow = c(G+1,1), mar = c(4,4,3,0.5))

   ylim <- range (sum_fit$scmuj)  
   for (g in 1:G)
   {
      plot (sum_fit$scmuj[,g], type = "h", ylim = ylim,
            ylab = "Normalized Mean",
            xlab = "Gene Rank by F-statistic",
            main = sprintf ("Normalized Means (Signals) of Class %d", g))
   }
   
   plot (sum_fit$signalj, type = "h", ylim = c(0, max(sum_fit$signalj)), 
         ylab = "Average Signal Level",
         xlab = "Gene Rank by F-statistic",
         main = "Overall Signal Levels of Top Genes")
}

############################################################################
######################### Utility Functions ################################
############################################################################

## compute log (sum (exp (lx) ))
log_sum_exp <- function(lx)
{  mlx <- max(lx)
   log(sum(exp(lx - mlx))) + mlx
}

## draw random numbers from inverse gamma distribution
rinvgam <- function (n, alpha, lambda)
{
    1 / rgamma (n, alpha, 1) * lambda
}

richisq <- function (n, alpha, w = 1)
{
    1 / rgamma (n, alpha / 2) * alpha * w / 2
}

## this is a generic function for generating Markov chain samples
## from a given density with Metropolis method
met_gauss <- function
   ( iters = 100, log_f, ini_value, stepsize = 0.5, diag_mh = FALSE, ...)
{
    state <- ini_value
    no_var <- length (state)
    mchain <- matrix (state, no_var , iters)
    nos_rej <- 0
    logf <- log_f (state,...)

    if (!is.finite (logf)) stop ("Initial value has 0 probability")

    for (i in 1:iters)
    {
        new_state <- rnorm (no_var, state, stepsize)
        new_logf <- log_f (new_state,...)

        if (log (runif(1)) < new_logf - logf)
        {
            state <- new_state
            logf <- new_logf
        }
        else nos_rej <- nos_rej + 1

        ## save state in chain
        mchain[,i] <- state
    }

    if (diag_mh)
    {
      cat ("Markov chain is saved in 'mchain' with columns for iterations\n")
      cat ("Rejection rate = ", nos_rej / iters, "\n")
      browser () ## pause to allow user to look at Markov chain
    }
    state
}


## This function estimates the parameters and hyerparameters based on
## the mle of means and variances of each feature.
trpr_mle <- function (X_tr, y_tr, X_ts = NULL,
            nos_fsel = ncol (X_tr), rankf = FALSE)
{
    ## read information about data
    n <- nrow (X_tr) ## numbers of obs
    p <- ncol (X_tr)
    ## find number of observations in each group
    nos_g <- as.vector (tapply (rep(1,n),INDEX = y_tr, sum))
    G <- length (nos_g)
    if (any(nos_g < 2)) stop ("Less than 2 cases in some group")
    freqy = nos_g / sum (nos_g)

    ## feature selection
    if (any (c(rankf, nos_fsel < p)) )
      info_sel <- rank_F (X_tr, y_tr)
    else
    { info_sel <- list (vars = 1:p)
    }

    ## create result storage
    nnfsel <- length (nos_fsel)

    list_fit_mle <- rep (list (""), nnfsel)

    array_probs_pred <- NULL
    if (!is.null (X_ts))
        array_probs_pred <- array (0, dim = c (nrow (X_ts), G, nnfsel) )

    for (i in 1:nnfsel)
    {
      k <- nos_fsel [i]

      if (k == 0)
      {
       fsel <- NULL
       fit_mle <- list (freqy = freqy)
      }
      else
      {

        fsel <- info_sel $ vars [1:k]

        X_tr_sel <- X_tr [, fsel, drop = FALSE]

        ## sufficient statistic and pooled variances
        gsum_X <- rowsum (X_tr_sel,y_tr)
        sum_X2 <- colSums (X_tr_sel^2)
        sum_gsum2 <- colSums (gsum_X^2 / nos_g )
        pvars <- (sum_X2 - sum_gsum2) / (n-G) ## pooled variances

        muj <- t(gsum_X / nos_g) ## group means
        wxj <- pvars ## pooled variances
        wx <- 1/mean (1/wxj)
        nuj <- rowMeans (muj)
        wnu <- mean (nuj^2)
        cmuj <- muj - nuj
        wmuj <- rowMeans (cmuj^2)
        wmu <- 1/mean (1/wmuj)
        scmuj <- cmuj/sqrt (wxj)

        fit_mle <- list (
            muj = muj, wxj = wxj, wmuj = wmuj, nuj = nuj, cmuj = cmuj,
            wx = wx, wmu = wmu, wnu = wnu, scmuj = scmuj,
            freqy = freqy,
            fsel = fsel )
      }

      list_fit_mle [[i]] <- fit_mle

      ## note: in "muj", the row is for features, the column is for groups
      if (!is.null (X_ts))
      {
        array_probs_pred [,,i] <- mlepred (X_ts = X_ts, fit_mle = fit_mle)
      }
    }

    list (
      array_probs_pred = array_probs_pred,
      nos_fsel = nos_fsel, list_fit_mle = list_fit_mle)
}

############################################################################
######################### Functions for Feature Selection ##################
############################################################################

## This function ranks all features in terms of F-statistic.
rank_F <- function(X, y)
{
  ## This function computes the values of F-statistic of all the features.
  comp_fstats <- function (X,y)
  {
    n <- length (y)
    nos_g <- as.vector (tapply (rep (1,n), INDEX = y, sum))
    G <- length (nos_g)

    gsum_X <- rowsum (X,y)
    sum_X2 <- colSums (X^2)
    sum_gsum2 <- colSums (gsum_X^2 / nos_g )

    pvars <- (sum_X2 - sum_gsum2) / (n-G) ## pooled variances

    sum_X <- colSums (X)
    gvars <- (sum_gsum2 - sum_X^2 / n) / (G-1) ## variances btw groups
    ## F-statistic
    gvars / pvars
  }

  fstats <- comp_fstats (X, y)
  vars <- order (fstats, decreasing = TRUE)
  list (vars=vars, fstats = fstats [vars])
}

############################################################################
######################### adjustment factor ################################
############################################################################

## This function generates random  part of lambda of poisson
## distribution and values of CDF of central F distribution, which are
## needed in approximating adjustment factor --- the probability that the
## F-statistic of a feature is smaller than a threshold
gen_qflmd <-  function (y_tr, cut_F, alpha1_mu = 1, alpha1_x = 10,
              cut_qf = exp (-10), nos_sim = 1000)
{
  n <- length (y_tr)
  nos_g <- tapply (rep(1,n), INDEX = y_tr, sum)
  G <- length(nos_g)

  qf <- c()
  l <- 1
  while (TRUE)
  {
    qf [l] <- pf(cut_F*(G-1)/(G-1 + 2*(l-1)), G-1 + 2*(l-1), n-G)
    if( qf[l] <= cut_qf ) break
    l <- l + 1
  }

  gen_adev <- function ()
  {
      mu <- rnorm (G)
      mu.bar <- sum(mu * nos_g) / n
      sum (mu^2 * nos_g) - n * mu.bar^2
  }

  devs <- replicate (nos_sim, gen_adev() )
  plmd <- devs/2 * richisq (nos_sim,alpha1_mu) / richisq (nos_sim,alpha1_x)

  list (qf = qf, plmd = plmd)
}

## Given 'qf' and 'plmd' returned by 'gen_qflmd', this function computes
## the probability that the  F-statistic is smaller than a threshold
comp_adjfactor <- function(w_mu, w_x, qflmd, cut_dpoi = exp (-10) )
{
  lmd <- qflmd$plmd * w_mu / w_x
  qf <- qflmd$qf
  .C("comp_adjfactor", PACKAGE = "BCBCSF",
  cut_dpoi, length(qf),length(lmd), qf, lmd, adjf = 0.0 )$adjf
}


