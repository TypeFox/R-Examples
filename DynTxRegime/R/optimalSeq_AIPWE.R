#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# Objective function for AIPWE                                                 #
# called only by rgenoud method for sequential DTR method                      #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# eta     : array of current parameters estimates for tx regime.               #
#                                                                              #
# regimes : an object of class Regimes or RegimesList                          #
#           For each dp (dp), a function defining the tx                       #
#           rule. For example, if the tx rule is : I(eta_1 < x1),              #
#           regimes is defined as                                              #
#             regimes <- function(a,data){as.numeric(a < data$x1)}             #
#           THE LAST ARGUMENT IS ALWAYS TAKEN TO BE THE DATA.FRAME             #
#                                                                              #
# txInfo  : an object of class TxInfo or TxInfoList                            #
#                                                                              #
# l.data  : data.frame of covariates and tx histories                          #
#                                                                              #
# outcome : an object of class QLearnEst or QLearnEstList.                     #
#                                                                              #
# propen  : an object of class PropenFit or PropenFitList                      #
#                                                                              #
# moMain  : an object of class modelObj, a list of objects of class modelObj,  #
#           or a list of objects of class modelObjSubset, which define the     #
#           models and R methods to be used to obtain parameter estimates and  #
#           predictions for for the main effects component of the              #
#           outcome regression.                                                #
#           See ?modelObj and/or ?modelObjSubset for details.                  #
#                                                                              #
# moCont  : an object of class modelObj, a list of objects of class modelObj,  #
#           or a list of objects of class modelObjSubset, which define the     #
#           models and R methods to be used to obtain parameter estimates and  #
#           predictions for for the contrasts component of the                 #
#           outcome regression.                                                #
#           See ?modelObj and/or ?modelObjSubset for details.                  #
#                                                                              #
# refit   : TRUE/FALSE - refit conditional expectations at each new set of     #
#           tx regime parameters.                                              #
#                                                                              #
# response: Response vector                                                    #
#                                                                              #
# iter    : an integer                                                         #
#                                                                              #
#           >=1 if moMain and moCont are to be fitted iteratively              #
#           The value is the maximum number of iterations.                     #
#           Note the iterative algorithms is as follows:                       #
#           Y = Ymain + Ycont                                                  #
#            (1) hat(Ycont) = 0                                                #
#            (2) Ymain = Y - hat(Ycont)                                        #
#            (3) fit Ymain ~ moMain                                            #
#            (4) set Ycont = Y - hat(Ymain)                                    #
#            (5) fit Ycont ~ moCont                                            #
#            (6) Repeat steps (2) - (5) until convergence or                   #
#            a maximum of iter iterations.                                     #
#                                                                              #
#           <=0 moMain and moCont will be combined and fit as a single object. #
#                                                                              #
#           Note that if iter <= 0, all non-model components of the            #
#           moMain and moCont must be identical                                #
#                                                                              #
# base    : An integer indicating the base tx or NULL (ordinal tx)             #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
optimalSeq_AIPWE <- function(eta, 
                             regimes, 
                             txInfo, 
                             l.data,  
                             outcome,  
                             propen,
                             moMain,
                             moCont,
                             refit,
                             response,
                             iter){
  nSamples <- nrow(l.data)

  if( is(regimes,'RegimeList') ) {
    nDP <- length(regimes)
    if( !is(txInfo,'TxInfoList') || length(txInfo) != nDP) {
      DeveloperError("Length mismatch between regimes and txInfo.",
                     "optimalSeq_AIPWE")
    }
  } else if( is(regimes,'Regime') ){
    nDP <- 1L
    refit <- FALSE
    if( !is(txInfo,'TxInfo') ) {
      DeveloperError("Length mismatch between regimes and txInfo.",
                     "optimalSeq_AIPWE")
    }
  } else {
    DeveloperError(paste("Unknown class for regimes",
                         paste(is(regimes),collapse=","),sep=""),
                   "optimalSeq_AIPWE")
  }

  valueFunc <- matrix(data = 0.0, 
                      nrow = nSamples,  
                      ncol = nDP)

  #--------------------------------------------------------------------------#
  # ind indicates if patient tx is in accordance with regime at dp i.        #
  # lambda : Probability that the tx does not follow regime. Pr(A_k != g_i)  #
  # leta : list of eta values used to eliminate excess as.list calls         #
  #--------------------------------------------------------------------------#
  ind <- matrix(data = 0L,  
                nrow = nSamples,  
                ncol = nDP)
  lambda <- matrix(data = 0.0,  
                   nrow = nSamples,  
                   ncol = nDP)
  leta <- as.list(eta)

  #--------------------------------------------------------------------------#
  # DR is a vector of values for the doubly robust estimator of each patient #
  #--------------------------------------------------------------------------#
  DR <- array(data = 0.0, dim = nSamples)

  #--------------------------------------------------------------------------#
  # Cycle through each dp to obtain probabilities, indicators and fits       #
  # j : local variable points to next unused eta parameter in leta           #
  #--------------------------------------------------------------------------#
  j <- length(eta)

  for( i in nDP:1L ){
    #----------------------------------------------------------------------#
    # For the current tx rule parameter values, obtain tx for each patient #
    # per the rule.                                                        #
    # reg.g : tx per regime                                                #
    #----------------------------------------------------------------------#
    if( nDP == 1L ) {
      txNm <- TxName(txInfo)
      rgm <- regimes
      txObj <- txInfo
    } else {
      txNm <- TxName(txInfo[[i]])
      rgm <- regimes[[i]]
      txObj <- txInfo[[i]]
    }

    #----------------------------------------------------------------------#
    # Number of parameters in regime i                                     #
    #----------------------------------------------------------------------#
    nPar <- NVars(rgm)

    #----------------------------------------------------------------------#
    # Create argument list for regime function                             #
    #----------------------------------------------------------------------#
    argList <- c(leta[(j-nPar+1L):j])
    names(argList) <- VNames(rgm)[1L:nPar]
    argList[[ VNames(rgm)[-c(1L:nPar)] ]] <- quote(l.data)

    #----------------------------------------------------------------------#
    # Calculate treatment regime values at current parameter estimates     #
    #----------------------------------------------------------------------#
    reg.g <- do.call(what = RegFunc(rgm), args = argList)

    #----------------------------------------------------------------------#
    # Convert to appropriate class                                         #
    #----------------------------------------------------------------------#
    if( is(l.data[,txNm], "factor") ) {
      reg.g <- factor(reg.g, levels = levels(l.data[,txNm]))
    } else {
      reg.g <- as.integer(round(reg.g,0L))
    }

    #----------------------------------------------------------------------#
    # Verify that value returned by regime is allowed by fSet              #
    #----------------------------------------------------------------------#
    for( ti in 1L:length(txObj@subsets) ) {
      inss <- txObj@ptsSubset %in% names(txObj@subsets)[ti]
      tst <- reg.g[inss] %in% txObj@subsets[[ti]]
      if( !all(tst) ) {
        stop("regime returns a value not allowed by fset. Verify inputs.")
      }
    }

    #----------------------------------------------------------------------#
    # move j to point to next unused eta value in leta                     #
    #----------------------------------------------------------------------#
    j <- j - nPar

    #----------------------------------------------------------------------#
    # ind[,i] = 1 if patient tx in accordance with regime at dp i.         #
    #           I(A_i = g_i)                                               #
    #         = 0 if patient tx not in accordance with regime at dp i.     #
    #           I(A_i != g_i)                                              #
    #----------------------------------------------------------------------#
    ind[,i] <- as.integer(l.data[,txNm] == reg.g | 
                          (is.na(l.data[,txNm]) & is.na(reg.g)))

    #----------------------------------------------------------------------#
    # If requested, refit conditional expectation models at new set of     #
    # tx regime parameter values.                                          #
    #----------------------------------------------------------------------#
    if(refit){
      if( i != nDP ){
        res <- valueFunc[,{i+1L}]
      } else {
        res <- response
      }

      if( nDP > 1L ){
        txI <- txInfo[[i]]
      } else {
        txI <- txInfo
      }

      mods <- pickOutcomeModels(moMain, moCont, i)

      tempOutcome <- qLearnEst(moMain = mods$moMain,
                               moCont = mods$moCont,
                               data = l.data,
                               response = res,
                               iter = iter, 
                               txInfo = txI)

      #------------------------------------------------------------------#
      # Eliminate patients with only 1 tx option from dataset for fit    #
      #------------------------------------------------------------------#
      use4fit <- eliminateSingleTx(subsets=txI@subsets, 
                                   ptsSubset=txI@ptsSubset)

      #------------------------------------------------------------------#
      # For patients with only 1 tx, use response as the value function  #
      #------------------------------------------------------------------#
      valueFunc[,i] <- res

      #------------------------------------------------------------------#
      # Change ith tx for all patients to be in accordance with current  #
      # regime                                                           #
      #------------------------------------------------------------------#
      l.data[,txNm] <- reg.g

      #------------------------------------------------------------------#
      # Calculate value function                                         #
      #------------------------------------------------------------------#
      me <- PredictMain(tempOutcome, l.data[use4fit,,drop=FALSE])
      cn <- PredictCont(tempOutcome, l.data[use4fit,,drop=FALSE])

      valueFunc[use4fit,i]  <- me + cn

      YTilde(tempOutcome) <- valueFunc[,i]

      if( nDP > 1L ) {
        outcome[[i]] <- tempOutcome
      } else {
        outcome <- tempOutcome
      }

    } else {
      #------------------------------------------------------------------#
      # Change ith tx for all patients to be in accordance with current  #
      # regime                                                           #
      #------------------------------------------------------------------#
      l.data[,txNm] <- reg.g

    }

  }

  #--------------------------------------------------------------------------#
  # For each decision point...                                               #
  #--------------------------------------------------------------------------#
  for( i in 1L:nDP ){

    if( is(txInfo,'TxInfoList') ) {
      txI <- txInfo[[i]]
      proI <- propen[[i]]
    } else if( is(txInfo,'TxInfo') ) {
      txI <- txInfo
      proI <- propen
    } else {
      DeveloperError(paste("Unknown class for txInfo",
                           paste(is(txInfo),collapse=","),sep=""),
                     "optimalSeq_AIPWE")
    }

    #----------------------------------------------------------------------#
    # For patients with only one treatment option, do not predict          #
    #----------------------------------------------------------------------#
    use4fit <- eliminateSingleTx(subsets=txI@subsets, 
                                 ptsSubset=txI@ptsSubset)

    #----------------------------------------------------------------------#
    # Calculate propensity for treatment for subset of patients            #
    #----------------------------------------------------------------------#
    mm <- PredictPropen(object = proI, 
                        newdata = l.data[use4fit,,drop=FALSE])

    #----------------------------------------------------------------------#
    # Convert assigned treatment to character                              #
    #----------------------------------------------------------------------#
    if( is(l.data[,txI@txName],"factor") ) {
      reg.g <- levels(l.data[,txI@txName])[l.data[,txI@txName]]
      reg.g <- reg.g[use4fit]
    } else {
      reg.g <- as.character(round(l.data[use4fit,txI@txName],0L))
    }

    #----------------------------------------------------------------------#
    # Retrieve probability of patient being assigned treatment             #
    #----------------------------------------------------------------------#
    probOfG <- numeric(nrow(l.data))
    probOfG[use4fit] <- sapply(1L:sum(use4fit), function(x){mm[x,reg.g[x]]})

    #----------------------------------------------------------------------#
    # For those that have only 1 treatment, probability is 1.              #
    #----------------------------------------------------------------------#
    probOfG[!use4fit] <- 1.0

    lambda[,i] <- 1.0 - probOfG

  }

  #--------------------------------------------------------------------------#
  # doubly robust estimator.                                                 #
  #--------------------------------------------------------------------------#

  AC <- array(data = 1.0, dim = nSamples)
  pc <- array(data = 1.0, dim = nSamples)
  cumInd <- array(data = 1L, dim = nSamples)
  ind <- cbind(cumInd,ind)

  for(i in nDP:1L) {

    if( is(txInfo,'TxInfo') ) {
      txInfoTemp <- txInfo
      outc <- outcome
    } else if( is(txInfo,'TxInfoList') ){
      txInfoTemp <- txInfo[[i]]
      outc <- outcome[[i]]
    } else {
      DeveloperError(paste("Unknown class for txInfo",
                           paste(is(txInfo),collapse=","),sep=""),
                     "optimalSeq_AIPWE")
    }

    #----------------------------------------------------------------------#
    # Eliminate patients with only 1 tx option from dataset for fit        #
    #----------------------------------------------------------------------#
    use4fit <- eliminateSingleTx(subsets=Subsets(txInfoTemp), 
                                 ptsSubset=PtsSubset(txInfoTemp))

    #----------------------------------------------------------------------#
    # For patients with 1 tx, set to original Ytilde value                 #
    #----------------------------------------------------------------------#
    if( i == nDP) {
      valueFunc[,i] <- response
    } else {
      valueFunc[,i] <- valueFunc[,{i+1L}]
    }

    #----------------------------------------------------------------------#
    # Calculate value function using optimal tx set by regime              #
    #----------------------------------------------------------------------#
    me <- PredictMain(outc, l.data[use4fit,])
    cn <- PredictCont(outc, l.data[use4fit,])
    valueFunc[use4fit,i]  <- me + cn

  }

  for(i in 1L:nDP) {

    #----------------------------------------------------------------------#
    # cumInd = 1 if patient followed tx regime up to the ith dp.           #
    #            I(C_{eta} >= i)                                           #
    #        = 0 if patient did not follow tx regime up to the ith dp.     #
    #            I(C_{eta} < i)                                            #
    #----------------------------------------------------------------------#
    cumInd <- cumInd*ind[,i]

    #----------------------------------------------------------------------#
    # Ultimately,                                                          #
    # AC = 1 if all txs given to a patient follow the current regime.      #
    #        I(C_{eta} = infinity)                                         #
    #    = 0 otherwise.                                                    #
    #        I(C_{eta} <= K)                                               #
    #----------------------------------------------------------------------#
    AC <- AC*as.numeric(ind[,{i+1L}])

    #----------------------------------------------------------------------#
    # C = 1 if patient treated in accordance with regime up to dp i, but   #
    #       did not follow tx regime at dp i. I(C_{eta} = i)               #
    #   = 0 otherwise. I(C_{eta} != i)                                     #
    #----------------------------------------------------------------------#
    C <- cumInd*(1L - ind[,{i+1L}])

    #----------------------------------------------------------------------#
    # pc = probability that coarsening occurs at a later dp.               #
    #      Pr(C_{et} > i) = prod_{k=1}^{i} (Pr(A_k=g_k))                   #
    #                       prod_{k=1}^{i} (1-Pr(A_k!=g_k))                #
    #----------------------------------------------------------------------#
    pc <- pc*(1.0 - lambda[,i])

    #----------------------------------------------------------------------#
    #          I(C_{eta} = i) - Pr(A_i != g_i)*I(C_{eta} >= i)             #
    # DR = sum -----------------------------------------------  mu_i       #
    #       i               Pr(C_{eta} > i)                                #
    #----------------------------------------------------------------------#
    DR <- DR + (as.numeric(C) - lambda[,i]*as.numeric(cumInd))/pc *
                as.vector(valueFunc[,i])
  }

  #--------------------------------------------------------------------------#
  #     (   I(C_{eta} = infinity)        )                                   #
  # mean|   --------------------- Y + DR |                                   #
  #     (      Pr(C_{eta} > K)           )                                   #
  #--------------------------------------------------------------------------#
  mn <- sum(DR + AC/pc*as.vector(response), na.rm=TRUE)/as.numeric(nSamples)

  return(mn)
}
