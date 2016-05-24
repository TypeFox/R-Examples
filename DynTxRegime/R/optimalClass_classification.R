#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# optimalClass_classification : prepares data and calls classification method  #
# specified by user.                                                           #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# contrast : Vector of the value of the contrast function for each sample.     #
#                                                                              #
# moClass  : an object of class modelObj, which defines the models and         #
#            R methods to be used to obtain parameter estimates and            #
#            predictions for the classification                                #
#                                                                              #
#            It is assumed that the solver.method contains                     #
#              weights : A vector of weights to be used in the fitting         #
#                        process.                                              #
#                                                                              #
# data     : data frame of covariates and response                             #
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= Returns a list                                                             =#
#=   cf: classification fit object                                            =#
#=  opt: optimal treatment regime for training set                            =#
#=                                                                            =#
#==============================================================================#
optimalClass_classification <- function(contrast, 
                                        moClass,  
                                        data){

  #--------------------------------------------------------------------------#
  # Classification weight variable.                                          #
  #--------------------------------------------------------------------------#
  weights <- abs(contrast)

  #--------------------------------------------------------------------------#
  # Normalize weights                                                        #
  #--------------------------------------------------------------------------#
  norm.weights <- weights/sum(weights)

  #--------------------------------------------------------------------------#
  # Add weights to formal arguments of classification method                 #
  #--------------------------------------------------------------------------#
  cArgs <- solverArgs(moClass)
  cArgs[[ "weights" ]] <- norm.weights
  solverArgs(moClass) <- cArgs

  #--------------------------------------------------------------------------#
  # Classification labels                                                    #
  #--------------------------------------------------------------------------#
  ZinternalZ <- as.numeric(contrast > -1.5e-8)
  ZinternalZ <- as.factor(ZinternalZ)

  #--------------------------------------------------------------------------#
  # Obtain classification fit using fit routine of modelObj                  #
  #--------------------------------------------------------------------------#
  cf <- fit(object = moClass, 
            data = data, 
            response = ZinternalZ)

  #--------------------------------------------------------------------------#
  # Predict classification of training set                                   #
  #--------------------------------------------------------------------------#
  pred <- factor(predict(object = cf, newdata = data), levels=c("0","1"))

  return(list("cf" = cf, "opt" = pred))
}

