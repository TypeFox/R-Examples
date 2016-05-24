#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# optimalSeq_IPWE : Objective function for IPWE                                #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# eta      : array of current parameters estimates for tx regime.              #
#                                                                              #
# regimes  : a function or a list of functions.                                #
#            For each dp (dp), a function defining the tx                      #
#            rule. For example, if the tx rule is : I(eta_1 < x1),             #
#            regimes is defined as                                             #
#              regimes <- function(a,data){as.numeric(a < data$x1)}            #
#            THE LAST ARGUMENT IS ALWAYS TAKEN TO BE THE DATA.FRAME            #
#                                                                              #
# txInfo   : an object of class txInfo or a list of objects of class txInfo    #
#            Each element of the list corresponds to a single dp               #
#              @column : column index of data that contains the tx variable    #
#              @name : the tx variable column header in data                   #
#              @options is a vector of all tx options available at the dp.     #
#                                                                              #
# l.data   : data.frame of covariates                                          #
#                                                                              #
# propen   : object of class PropenFit or PropenFitList                        #
#                                                                              #
# response : Response vector                                                   #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
optimalSeq_IPWE <- function(eta, 
                            regimes,  
                            txInfo,  
                            l.data,  
                            propen,
                            response){

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

    } else {

      txNm <- TxName(txInfo[[i]])
      rgm <- regimes[[i]]

    }

    if( is(l.data[,txNm], "factor") ) {
      tdata <- l.data
      tdata[,txNm] <- levels(tdata[,txNm])[tdata[,txNm]]
    } else {
      tdata <- l.data
    }

    nPar <- NVars(rgm)
    argList <- c(leta[(j-nPar+1L):j])
    names(argList) <- VNames(rgm)[1L:nPar]
    argList[[ VNames(rgm)[-c(1L:nPar)] ]] <- quote(tdata)

    reg.g <- do.call(what = RegFunc(rgm), args = argList)

    if( is(l.data[,txNm], "factor") ) {
      reg.g <- factor(reg.g, levels = levels(l.data[,txNm]))
    } else {
      reg.g <- as.integer(round(reg.g,0L))
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
    # Change ith tx for all patients to be in accordance with current      #
    # regime                                                               #
    #----------------------------------------------------------------------#
    l.data[,txNm] <- reg.g


  }

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
    useGrps <- sapply(X = txI@subsets, FUN = length) > 1.5
    use4fit <- txI@ptsSubset %in% names(txI@subsets)[which(useGrps)]

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
    # Retrieve probability of treatment being assigned treatment           #
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
  # estimator.                                                               #
  #--------------------------------------------------------------------------#

  AC <- array(data = 1.0, dim = nSamples)
  pc <- array(data = 1.0, dim = nSamples)
  cumInd <- array(data = 1L, dim = nSamples)
  ind <- cbind(cumInd,ind)


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

  }

  #--------------------------------------------------------------------------#
  #     (   I(C_{eta} = infinity)        )                                   #
  # mean|   --------------------- Y + DR |                                   #
  #     (      Pr(C_{eta} > K)           )                                   #
  #--------------------------------------------------------------------------#
  mn <- sum(AC/pc*as.vector(response), na.rm=TRUE)/as.numeric(nSamples)

  return(mn)
}
