#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# optimalSeq : Main function call for robust estimation of optimal dynamic     #
#  tx regimes for sequential tx decisions as described in                      #
#                                                                              #
#    Baqun Zhang, Anastasios A. Tsiatis, Eric B. Laber & Marie Davidian,       #
#    "A Robust Method for Estimating Optimal Treatment Regimes", Biometrics    #
#    68, 1010-1018.                                                            #
#                                                                              #
#    and                                                                       #
#                                                                              #
#    Baqun Zhang, Anastasios A. Tsiatis, Eric B. Laber & Marie Davidian,       #
#    "Robust estimation of optimal treatment regimes for sequential treatment  #
#    decisions", Biometrika (2013) pp.1-14.                                    #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# moPropen: an object of class modelObj, a list of objects of class modelObj,  #
#           or a list of objects of class modelObjSubset, which define the     #
#           models and R methods to be used to obtain parameter estimates and  #
#           predictions for the propensity for treatment.                      #
#                                                                              #
#           If the prediction method specified in moPropen returns             #
#           predictions for only a subset of the categorical tx data,          #
#           it is assumed that the base level defined by levels(tx)[1] is      #
#           the missing category.                                              #
#                                                                              #
#           Note that it is assumed that the columns of the predictions are    #
#           ordered in accordance with the vector returned by levels().        #
#                                                                              #
# moMain  : an object of class modelObj, a list of objects of class modelObj,  #
#           or a list of objects of class modelObjSubset, which define the     #
#           models and R methods to be used to obtain parameter estimates and  #
#           predictions for the main effects component of the outcome          #
#           regression.                                                        #
#                                                                              #
# moCont  : an object of class modelObj, a list of objects of class modelObj,  #
#           or a list of objects of class modelObjSubset, which define the     #
#           models and R methods to be used to obtain parameter estimates and  #
#           predictions for the contrasts component of the outcome             #
#           regression.                                                        #
#                                                                              #
# data    : a data frame of the covariates and tx histories                    #
#           tx variableS will be recast as factors if not provided as such.    #
#                                                                              #
# response: response vector                                                    #
#                                                                              #
# txName  : a vector of characters.                                            #
#           The column headers of \emph{data} that correspond to the tx        #
#           covariate for each decision point.                                 #
#           The ordering should be sequential, i.e., the 1st element           #
#           gives column name for the 1st decision point tx, the               #
#           2nd gives column name for the 2nd decision point tx,               #
#           etc.                                                               #
#                                                                              #
# regimes : a function or a list of functions.                                 #
#           For each decision point, a function defining the tx                #
#           rule. For example, if the tx rule is : I(eta_1 < x1),              #
#           regimes is defined as                                              #
#             regimes <- function(a,data){as.numeric(a < data$x1)}             #
#           THE LAST ARGUMENT IS ALWAYS TAKEN TO BE THE DATA.FRAME             #
#                                                                              #
# fSet    : A function or a list of functions.                                 #
#           This argument allows the user to specify the subset of tx          #
#           options available to a patient or the subset of patients that      #
#           will be modeled uniquely.                                          #
#           The functions should accept as input either                        #
#             1) explicit covariate names as given in column headers of data   #
#             2) a vector of covariates (i.e. a row of a data.frame)           #
#           and must return a list. The first element of the list is a         #
#           character giving a nickname to the subset. The second element      #
#           is a vector of the tx options available to the subset.             #
#                                                                              #
# refit   : TRUE/FALSE flag indicating if outcome regression models            #
#           are to be refit for each new set of tx regime                      #
#           parameters (TRUE), or if Q-learning is to be used (FALSE).         #
#           Can only be 'true' for multiple decision point analyses.           #
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
# ...     : Additional arguments required by rgenoud. At a minimum this        #
#           should include Domains, pop.size and starting.values.              #
#           See ?rgenoud for more information.                                 #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
optimalSeq <- function(...,
                       moPropen,
                       moMain,
                       moCont,
                       data,
                       response,
                       txName,
                       regimes,
                       fSet = NULL,
                       refit = FALSE,
                       iter = 0L,
                       suppress = FALSE){

  if( !requireNamespace("rgenoud", quietly=TRUE) ) {
    stop("optimalSeq requires the rgenoud package available on CRAN.")
  }

  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  #++++++                         Verify Input                         ++++++#
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  #--------------------------------------------------------------------------#
  # Store call for return list                                               #
  #--------------------------------------------------------------------------#
  call <- match.call()

  #--------------------------------------------------------------------------#
  # Convert regimes input into local class definition regime or regimeList   #
  #--------------------------------------------------------------------------#
  if( !is.list(regimes) ){
    if( is.function(regimes) ){
      nDP <- 1L
    } else {
      UserError("input",
                "regime must be provided as a function.")
    }

    #----------------------------------------------------------------------#
    # Extract the formal arguments of the user function.                   #
    #----------------------------------------------------------------------#
    nms <- names(formals(regimes))
    #----------------------------------------------------------------------#
    # Count the number of argument to user function.                       #
    #----------------------------------------------------------------------#
    nVars <- as.integer(round(length(nms),0L)) - 1L

    #----------------------------------------------------------------------#
    # If function has no arguments or data is not in the list of names,    #
    # throw error.                                                         #
    #----------------------------------------------------------------------#
    if( nVars <= 0L || !('data' %in% nms) ) {
      UserError("input",
                paste("Formal arguments of function input through regimes ",
                 "must contain all regime parameters and 'data'.", sep=""))
    }

    #----------------------------------------------------------------------#
    # Create object of class Regime.                                       #
    #----------------------------------------------------------------------#
    regimes <- new("Regime",
                   nVars = as.integer(nVars),
                   vNames = as.vector(nms),
                   func = regimes )
  } else {
    nDP <- as.integer(round(length(regimes),0L))

    regimeInfo <- list()
    nVars <- 0L

    for( i in 1L:nDP ){

      if( !is.function(regimes[[i]]) ) {
        UserError("input",
                  "Each element of regimes must be a function.")
      }

      #------------------------------------------------------------------#
      # Extract the formal arguments of the user function.               #
      #------------------------------------------------------------------#
      nms <- names(formals(regimes[[i]]))

      #------------------------------------------------------------------#
      # Count the number of argument to user function.                   #
      #------------------------------------------------------------------#
      tvars <- as.integer(round(length(nms),0L)) - 1L
      #------------------------------------------------------------------#
      # If function has no arguments or data is not in the list of names,#
      # throw error.                                                     #
      #------------------------------------------------------------------#
      if( tvars <= 0L || !('data' %in% nms) ) {
        UserError("input",
                  paste("Formal arguments of function input through regimes ",
                        "must contain all regime parameters and 'data'.", 
                        sep=""))
      }

      #------------------------------------------------------------------#
      # Track total number of variables in all regimes.                  #
      #------------------------------------------------------------------#
      nVars <- nVars + tvars

      #------------------------------------------------------------------#
      # Create object of class Regime.                                   #
      #------------------------------------------------------------------#
      regimeInfo[[i]] <- new("Regime",
                             nVars = as.integer(tvars),
                             vNames = as.vector(nms),
                             func = regimes[[i]] )
    }

    #----------------------------------------------------------------------#
    # Create object of class RegimeList.                                   #
    #----------------------------------------------------------------------#
    regimes <- new("RegimeList", 
                   loo = regimeInfo)
    rm(regimeInfo)
  }


  #--------------------------------------------------------------------------#
  # Identify the type of estimator requested.                                #
  #--------------------------------------------------------------------------#
  if( is(moCont, "NULL") && is(moMain, "NULL") ) {
    if( !suppress ) cat("Inverse Probability Weighted Estimator\n")
    method <- 'ipwe'
  } else {
    if( !suppress ) cat("Augmented Inverse Probability Weighted Estimator\n")
    method <- 'aipwe'
  }

  if( is(moMain, "list") && all(is(unlist(moMain),"NULL")) ) {
    moMain <- NULL
  }

  if( is(moCont, "list") && all(is(unlist(moCont),"NULL")) ) {
    moCont <- NULL
  }

  #--------------------------------------------------------------------------#
  # For each modeling object, verify the number of models provided.          #
  # Convert to lists if multi-decision point.                                #
  #--------------------------------------------------------------------------#
  objs <- list("moPropen" = moPropen,
               "moMain" = moMain,
               "moCont" = moCont)

  for( i in 1L:3L ) {
    if( is(objs[[i]], "NULL") ) next

    if( is(objs[[i]],'modelObj') ){

      if(nDP != 1L) {
        UserError("input",
                  paste("Insufficient number of models specified for ",
                        names(objs)[i], ".", sep=""))
      }

    } else if( is(objs[[i]], "list") && !is(objs[[i]][[1]], "NULL") ) {

      if( is(objs[[i]][[1L]], 'modelObj') ) {
        clx <- 'ModelObj'
        cll <- 'modelObj'
      } else if( is(objs[[i]][[1L]],'ModelObjSubset') ) {
        clx <- 'ModelObjSubset'
        cll <- 'ModelObjSubset'
      }

      tst <- sapply(X = objs[[i]], FUN = class) != cll
      if( any(tst) ) {
        UserError("input",
                  paste("All elements of ", names(objs)[i],
                        " must be of the same class.", sep=""))
      }

      if(length(tst) < nDP) {
        UserError("input",
                  paste("Insufficient number of models specified for ",
                        names(objs)[i], ".", sep=""))
      }
      clxname <- paste(clx,"List",sep="")

      objs[[i]] <- new(clxname,
                       loo = objs[[i]])
    } else if( !is(objs[[i]][[1]], "NULL") ) {

      stop(paste(names(objs)[i], " must be an object of class modelObj, ",
                 "a list of objects of class modelObj, or ",
                 "a list of objects of class modelObjSubset.", sep=""))
    }
  }

  moPropen <- objs[[1L]]
  moMain <- objs[[2L]]
  moCont <- objs[[3L]]

  #--------------------------------------------------------------------------#
  # txName must be an object of class character                              #
  #--------------------------------------------------------------------------#
  if( !is(txName, "character") ) {
    UserError("input",
              "'txName' must be a character.")
  }

  #--------------------------------------------------------------------------#
  # Verify that the treatments are in dataset.                               #
  #--------------------------------------------------------------------------#
  for( i in 1L:length(txName) ) {

    txVec <- try(data[,txName[i]], silent = TRUE)

    if( is(txVec,"try-error") ) {
      UserError("input",
                paste(txName[i], " not found in data.", sep="") )
    }

    #----------------------------------------------------------------------#
    # If not a factor variable, convert to integer.                        #
    #----------------------------------------------------------------------#
    if( !is(txVec,"factor") ) {
      if( !isTRUE(all.equal(txVec, round(txVec,0L))) ) {
        UserError("input",
                  "Treatment variable must be a factor or an integer.")
      }
      data[,txName[i]] <- as.integer(round(data[,txName[i]],0L))
    }

  }

  response <- matrix(data = response,
                     ncol = 1L, 
                     dimnames = list(NULL,"YinternalY"))

  #--------------------------------------------------------------------------#
  # Process treatment and feasibility input                                  #
  #--------------------------------------------------------------------------#
  if( is(fSet, "list") ) {

    if( as.integer(round(length(fSet),0L)) != nDP ){
      UserError("input",
                "There must be one feasibility rule for each decision point.")
    }

    tst <- sapply(X=fSet, FUN=class) != 'function'
    if( any(tst) ) {
      UserError("input",
                "fSet must be a list of functions.")
    }

    txInfo <- list()
    for( i in 1L:nDP ){
      txInfo[[i]] <- txProcess(txVar = txName[i], 
                               data = data, 
                               fSet = fSet[[i]])
    }

    txInfo <- new("TxInfoList",
                  loo = txInfo)

  } else if( is(fSet, "function") ){

    if( nDP != 1L ) {
      UserError("input",
                "There must be one feasibility rule for each decision point.")
    }

    txInfo <- txProcess(txVar = txName, 
                        data = data, 
                        fSet = fSet)

  } else if( is(fSet,"NULL") ) {

    if( nDP == 1L ) {
      txInfo <- txProcess(txVar = txName, 
                          data = data, 
                          fSet = fSet)
    } else {
      txInfo <- list()
      for( i in 1L:nDP ){
        txInfo[[i]] <- txProcess(txVar = txName[i], 
                                 data = data, 
                                 fSet = fSet)
      }

      txInfo <- new("TxInfoList",
                    loo = txInfo)
    }

  } else {

    UserError("input",
              "If provided, fSet must be a function or a list of functions.")

  }

  #--------------------------------------------------------------------------#
  # Obtain fits for propensity and outcome                                   #
  # optimalCore returns a list.                                              #
  #  $outcome : a list of qLearnEst objects. If refitting is requested, a    #
  #             single object is the nDP slot                                #
  #  $propensity: a list of propenObj objects                                #
  #--------------------------------------------------------------------------#
  core <- optimalCore(moPropen = moPropen, 
                      moMain = moMain, 
                      moCont = moCont,
                      data = data, 
                      response = response, 
                      txInfo = txInfo, 
                      refit = refit, 
                      iter = iter)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#      Obtain genetic algorithm estimates for tx regime parameters             #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #----------------------------------------------------------------------#
  # create argument list specific to estimator requested by user.        #
  #----------------------------------------------------------------------#
  argList <- list(...)
  argList[['nvars']] <- nVars
  argList[['regimes']] <- regimes
  argList[['txInfo']] <- txInfo
  argList[['l.data']] <- quote(data)
  argList[['propen']] <- core$propen


  if(tolower(method)=='aipwe'){
    argList[['fn']] <- optimalSeq_AIPWE
    argList[['outcome']] <- core$outcome
    argList[['moMain']] <- moMain
    argList[['moCont']] <- moCont
    argList[['refit']] <- refit
    argList[['response']] <- response
    argList[['iter']] <- iter
  } else if(tolower(method)=='ipwe') {
    argList[['fn']] <- optimalSeq_IPWE
    argList[['response']] <- response
  }

  #--------------------------------------------------------------------------#
  # set argument list for rgenoud method                                     #
  #--------------------------------------------------------------------------#
  argList[['print.level']] <- 0L
  argList[['max']] <- TRUE
  argList[['gradient.check']] <- FALSE
  argList[['BFGS']] <- FALSE
  argList[['P9']] <- 0L
  argList[['optim.method']] <- "Nelder-Mead"

  #--------------------------------------------------------------------------#
  # Initiate genetic algorithm                                               #
  #--------------------------------------------------------------------------#
  gaEst <- do.call(rgenoud::genoud, argList)

  leta <- as.list(gaEst$par)
  j <- length(leta)

  optTx <- matrix(data = NA,
                  nrow = nrow(data),
                  ncol = nDP,
                  dimnames = list(NULL,paste("dp=", 1L:nDP)))

  #--------------------------------------------------------------------------#
  # Calculate optimal treatment regime and refit models if requested.        #
  #--------------------------------------------------------------------------#
  k <- nDP
  while( k > 0L ){

    if( is(regimes,'RegimeList') ) {
      rgm <- regimes[[k]]
    } else if( is(regimes,'Regime') ) {
      rgm <- regimes
    } else {
      DeveloperError(paste("regimes of wrong class", 
                           paste(is(regimes),collapse=","), sep=""),
                     "optimalSeq")
    }

    #----------------------------------------------------------------------#
    # For the current tx rule parameter values, obtain tx                  #
    # for each patient per the rule defined at the kth dp.                 #
    #----------------------------------------------------------------------#
    if( is(data[,txName[k]], "factor") ) {
      tdata <- data
      tdata[,txName[k]] <- levels(tdata[,txName[k]])[tdata[,txName[k]]]
    } else {
      tdata <- data
    }
    argList <- c(leta[(j-NVars(rgm)+1L):j])
    names(argList) <- VNames(rgm)[1L:NVars(rgm)]
    argList[[ 'data' ]] <- quote(tdata)
    temp <- do.call(RegFunc(rgm), argList)

    if( is(data[,txName[k]], "factor") ) {
      optTx[,k] <- factor(temp,levels=levels(data[,txName[k]]))
    } else {
      optTx[,k] <- as.integer(round(temp,0L))
    }

    #----------------------------------------------------------------------#
    # move j to point to next unused eta value in leta                     #
    #----------------------------------------------------------------------#
    j <- j - NVars(rgm)

    if(refit){
      #------------------------------------------------------------------#
      # Refit conditional expectation values at final tx regime values.  #
      #------------------------------------------------------------------#
      mods <- pickOutcomeModels(moMain, moCont, k)

      if( k == nDP ) {
        resp <- response
      } else {
        resp <- YTilde(core$outcome[[k+1L]])
      }

      if( nDP == 1L ) {
        txI <- txInfo
      } else {
        txI <- txInfo[[k]]
      }

      tempOutcome <- qLearnEst(moMain = mods$moMain,
                               moCont = mods$moCont,
                               data = data,
                               txName = TxName(txI),
                               response = resp,
                               iter = iter,
                               txInfo = txI)

      #------------------------------------------------------------------#
      # Eliminate patients with only 1 tx option from dataset for fit    #
      #------------------------------------------------------------------#
      use4fit <- eliminateSingleTx(subsets=Subsets(txI), 
                                   ptsSubset=PtsSubset(txI))

      #------------------------------------------------------------------#
      # Change kth tx for all patients to be in accordance with regime   #
      #------------------------------------------------------------------#
      data[,TxName(txI)] <- optTx[,k]

      #------------------------------------------------------------------#
      # Calculate value function                                         #
      #------------------------------------------------------------------#
      me <- PredictMain(tempOutcome, data[use4fit,,drop=FALSE])
      cn <- PredictCont(tempOutcome, data[use4fit,,drop=FALSE])

      yt <- YTilde(tempOutcome)
      yt[use4fit] <- me + cn
      YTilde(tempOutcome) <- yt
      if( nDP == 1L ) {
        core$outcome <- tempOutcome
      } else {
        core$outcome[[k]] <- tempOutcome
      }
    }
    k <- k - 1L
  }

  #--------------------------------------------------------------------------#
  # For each parameter of regime, add information to name indicating the dp  #
  #--------------------------------------------------------------------------#
  pars <- list()
  if( is(regimes,'RegimeList') ) {
    k <- 1L
    for(i in 1L:nDP){
      pars[[i]] <- c(gaEst$par[k:{k+NVars(regimes[[i]])-1L}])
      k <- k + NVars(regimes[[i]])
    }
  } else if( is(regimes,'Regime') ) {
    pars <- list(gaEst$par)
  } else {
      DeveloperError(paste("regimes of wrong class", 
                           paste(is(regimes),collapse=","), sep=""),
                     "optimalSeq")
  }


  #--------------------------------------------------------------------------#
  # Create new object of class optimal to be returned to user.               #
  #--------------------------------------------------------------------------#
  result <- new("OptimalSeq", 
                varEst = pars,
                estVal = gaEst$value,
                genetic = gaEst,
                propen = core$propensity,
                outcome = core$outcome,
                regimes = regimes,
                optTx = optTx,
                txInfo = txInfo,
                call = call)

  if( !suppress ) show(result)

  return(result)
}
