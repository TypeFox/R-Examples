#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# iqEst - Fit main effect and contrasts components of first and second stage   #
#         IQ-Learning algorithm.                                               #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# moMain   : an object of class modelObj that defines the models and R methods #
#            to be used to obtain parameter estimates and predictions for main #
#            effects component of the outcome regression.                      #
#                                                                              #
# moCont   : an object of class modelObj that defines the models and R methods #
#            to be used to obtain parameter estimates and predictions for the  #
#            contrasts component of the outcome regression.                    #
#                                                                              #
# data     : data.frame of covariates and treatment histories                  #
#                                                                              #
# response : response vector                                                   #
#                                                                              #
# txInfo   : object of class txInfo                                            #
#                                                                              #
# iter     : an integer                                                        #
#                                                                              #
#           >=1 if moMain and moCont are to be fit using iterative algorithm   #
#           The value of iter is the maximum number of iterations.             #
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
#           Either categorical or integer data can be provided for the tx.     #
#           If categorical, the fitted contrast and main effects are defined   #
#           relative to the base category {defined as levels()[1]}. The values #
#           may not be those returned by predict(object) if iterate fits are   #
#           used. If integer, the fitted contrast and main effects are defined #
#           relative to no tx (tx = 0).                                        #
#                                                                              #
#           Note that if iter <= 0, all non-model components of the            #
#           moMain and moCont must be identical                                #
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= Returns an object of class iqEst                                           =#
#=                                                                            =#
#==============================================================================#
iqEst <- function(moMain,
                  moCont,
                  data, 
                  response, 
                  txInfo,
                  iter, ...){


  #--------------------------------------------------------------------------#
  # fit model using provided solver                                          #
  #--------------------------------------------------------------------------#
  fitObj <- fitChoose(moMain = moMain, 
                      moCont = moCont, 
                      data = data,
                      response = response,
                      txName = TxName(txInfo), 
                      iter = iter)

  result <- new("IQEst",
                fitObj = fitObj)

  return(result)
}
