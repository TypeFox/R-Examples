#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# optimalClass : implementation of estimation of optimal dynamic tx            #
#                 regimes from a classification perspective.                   #
#                                                                              #
#    Baqun Zhang, Anastasios A. Tsiatis, Marie Davidian, Min Zhang and         #
#    Eric B. Laber. "Estimating optimal tx regimes from a classification       #
#    perspective." Stat 2012; 1: 103-114.                                      #
#                                                                              #
# Note that this method is a single decision point, binary treatment method    #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# moPropen: an object of class modelObj, which defines the models and R        #
#           methods to be used to obtain parameter estimates and               #
#           predictions for the propensity for treatment.                      #
#                                                                              #
#           If the prediction method specified in moPropen returns             #
#           predictions for only a subset of the categorical tx data,          #
#           it is assumed that the base level defined by levels(tx)[1] is      #
#           the missing category.                                              #
#                                                                              #
# moMain  : an object of class modelObj, which defines the models and R        #
#           methods to be used to obtain parameter estimates and               #
#           predictions for for the main effects component of the              #
#           outcome regression.                                                #
#           NULL is an appropriate value.                                      #
#                                                                              #
# moCont  : an object of class modelObj, which defines the models and R        #
#           methods to be used to obtain parameter estimates and               #
#           predictions for for the contrasts component of the                 #
#           outcome regression.                                                #
#           NULL is an appropriate value.                                      #
#                                                                              #
# moClass : an object of class modelObj, which defines the                     #
#           models and R methods to be used to obtain parameter                #
#           estimates and predictions for the classification.                  #
#                                                                              #
# data    : a data frame of the covariates and tx histories                    #
#           tx variable will be recast as factor if not provided as such.      #
#                                                                              #
# response: response vector                                                    #
#                                                                              #
# txName  : an character giving the column header of the column in data        #
#           that contains the tx covariate.                                    #
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
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= Returns an object of class optimalClass                                    =#
#=                                                                            =#
#==============================================================================#
optimalClass <- function(..., 
                         moPropen,
                         moMain,
                         moCont,
                         moClass,
                         data,
                         response,
                         txName,
                         iter = 0L,
                         suppress = FALSE){

  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  #++++++                         Verify Input                         ++++++#
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  #--------------------------------------------------------------------------#
  # Store call for return list                                               #
  #--------------------------------------------------------------------------#
  call <- match.call()

  #--------------------------------------------------------------------------#
  # iter must be an integer                                                  #
  #--------------------------------------------------------------------------#
  if( !is(iter, "integer") ) iter <- as.integer(round(iter,0L))

  #--------------------------------------------------------------------------#
  # If provided, modeling object must be single formula versions             #
  #--------------------------------------------------------------------------#
  if( !is(moMain, 'modelObj') && !is(moMain, 'NULL') ){
    UserError("input", 
              "moMain must be an object of class modelObj or NULL")
  }

  if( !is(moCont, 'modelObj') && !is(moCont, 'NULL') ){
    UserError("input", 
              "moCont must be an object of class modelObj or NULL")
  }

  #--------------------------------------------------------------------------#
  # Determine requested estimator.                                           #
  #--------------------------------------------------------------------------#
  if( is.null(moMain) && is.null(moCont) ){
    if( !suppress ) cat("\nInverse Probability Weighted Estimator\n")
    method <- "ipwe"
  } else {
    if( !suppress ) cat("\nAugmented Inverse Probability Weighted Estimator\n")
    method <- "aipwe"
  }

  #--------------------------------------------------------------------------#
  # moPropen must be provided and must be of class modelObj.                 #
  #--------------------------------------------------------------------------#
  if( !is(moPropen,"modelObj") ){
    UserError("input", 
              "moPropen must an object of class modelObj")
  }

  #--------------------------------------------------------------------------#
  # moClass must be provided and must be of class modelObj.                  #
  #--------------------------------------------------------------------------#
  if( !is(moClass,"modelObj") ){
    UserError("input", 
              "moClass must an object of class modelObj")
  }

  #--------------------------------------------------------------------------#
  # txName must be an object of class character                              #
  #--------------------------------------------------------------------------#
  if( !is(txName, "character") ) {
    UserError("input",
              "'txName' must be a character.")
  }

  #--------------------------------------------------------------------------#
  # verify that treatment variable is in provided data.frame.                #
  #--------------------------------------------------------------------------#
  txVec <- try(data[,txName], silent = TRUE)

  if( is(txVec,"try-error") ) {
    UserError("input",
              paste(txName, " not found in data.", sep="") )
  }

  #--------------------------------------------------------------------------#
  # Treatment must be an integer.                                            #
  # If a factor, throw error. If numeric, recast as integer.                 #
  #--------------------------------------------------------------------------#
  if( is(txVec,"factor") ) {
      UserError("input",
                "Treatment variable must be an integer.")
  } else {
    if( !isTRUE(all.equal(txVec, round(txVec,0L))) ) {
      UserError("input",
                "Treatment variable must be an integer.")
    }
    data[,txName] <- as.integer(round(data[,txName],0L))
  }

  #--------------------------------------------------------------------------#
  # Cast response as a matrix                                                #
  #--------------------------------------------------------------------------#
  response <- matrix(data = response,
                     ncol = 1L, 
                     dimnames = list(NULL,"YinternalY"))

  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  #++++++                         Calculation                          ++++++#
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  #--------------------------------------------------------------------------#
  # Process treatment information.                                           #
  #--------------------------------------------------------------------------#
  txInfo <- txProcess(txVar = txName, 
                      data = data, 
                      fSet = NULL)

  #--------------------------------------------------------------------------#
  # Verify that treatment if binary.                                         #
  #--------------------------------------------------------------------------#
  superSet <- SuperSet(txInfo)
  if( as.integer(round(length(superSet),0L)) != 2L ) {
    UserError("input",
              "Only binary treatment options can be used in this method.")
  }

  #--------------------------------------------------------------------------#
  # Verify that treatment is {0,1}                                           #
  #--------------------------------------------------------------------------#
  if( any( !(superSet %in% c("0","1")) ) ){
    UserError("input",
              "Treatments must be coded as 0 and 1")
  }

  #--------------------------------------------------------------------------#
  # Obtain fits for propensity and outcome                                   #
  #--------------------------------------------------------------------------#
  core <- optimalCore(moPropen = moPropen, 
                      moMain = moMain, 
                      moCont = moCont,
                      data = data, 
                      response = response, 
                      txInfo = txInfo, 
                      refit = FALSE, 
                      iter = iter)

  #--------------------------------------------------------------------------#
  # Obtain propensity predictions.                                           #
  #--------------------------------------------------------------------------#
  propenMatrix <- PredictPropen(object = core$propen,
                                newdata = data, 
                                subset = txInfo@superSet)

  #--------------------------------------------------------------------------#
  # Calculate contrast                                                       #
  #--------------------------------------------------------------------------#
  if( tolower(method) == "aipwe" ) {

    contrast <- optimalClass_AIPWE(outcome = core$outcome, 
                                   txInfo = txInfo, 
                                   propensity = propenMatrix, 
                                   data = data,
                                   response = response)

  } else if( tolower(method) == "ipwe" ) {

    contrast <- optimalClass_IPWE(txInfo = txInfo, 
                                  propensity = propenMatrix, 
                                  data = data,
                                  response = response)

  }

  #--------------------------------------------------------------------------#
  # Obtain classification fit                                                #
  #--------------------------------------------------------------------------#
  classFit <- optimalClass_classification(contrast = contrast$contrast, 
                                          moClass = moClass,  
                                          data = data)

  #--------------------------------------------------------------------------#
  # Calculate estimator.                                                     #
  #--------------------------------------------------------------------------#
  grp1 <- levels(classFit$opt)[classFit$opt] > 0.5
  estResponse <- sum(contrast$contrast[grp1])/nrow(data) + contrast$mean.mu0

  #--------------------------------------------------------------------------#
  # Create OptimalClass object and return to user.                           #
  #--------------------------------------------------------------------------#
  oc1 <- new("OptimalClass",
             estVal = estResponse,
             classif = classFit$cf,
             propen = core$propen,
             outcome = core$outcome,
             optTx = classFit$opt,
             call = call)

  if( !suppress ) show(oc1)

  return(oc1)

}
