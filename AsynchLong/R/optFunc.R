#----------------------------------------------------------------------#
# optFunc : Function that R's optim is minimizing. Taken to be U^2     #
#----------------------------------------------------------------------#
#                                                                      #
# pars      : Current parameter estimates                              #
#                                                                      #
# data.y    : Matrix of responses. It is assumed that the first column #
#             contains integer patient IDs, the second column contains #
#             the time of measurement, and the third column contains   #
#             the value of the measurement.                            #
#                                                                      #
# data.x    : Matrix of covariates. The columns contain only the values#
#             of the covariates.                                       #
#                                                                      #
# kernel    : a character. One of "TD" = time dependent,               #
#             "TI" = time-invariant or "LV" = last value.              #
#             This specifies how the kernel is used                    #
#                                                                      #
# lType     : a character. One of "identity", "log", "logistic".       #
#             Specifies the form of the link function.                 #
#                                                                      #
# tt        : If provided, a vector of times at which to evaluate the  #
#             kernel                                                   #
#                                                                      #
#----------------------------------------------------------------------#
#                                                                      #
# Returns scalar value U^2                                             #
#                                                                      #
#----------------------------------------------------------------------#
optFunc <- function(pars,
                    data.y, 
                    data.x,
                    kernel, 
                    lType,
                    xIs,
                    yIs,
                    nPatients) {

  res <- cmp_uFunc(pars = pars,
                   data.y = data.y, 
                   data.x = data.x,
                   kernel = kernel, 
                   lType = lType,
                   xIs = xIs,
                   yIs = yIs,
                   nPatients = nPatients)

  return(res %*% res)


}

#----------------------------------------------------------------------#
# doptFunc : Derivative of function that R's optim is minimizing.      #
#----------------------------------------------------------------------#
#                                                                      #
# pars      : Current parameter estimates                              #
#                                                                      #
# data.y    : Matrix of responses. It is assumed that the first column #
#             contains integer patient IDs, the second column contains #
#             the time of measurement, and the third column contains   #
#             the value of the measurement.                            #
#                                                                      #
# data.x    : Matrix of covariates. The columns contain only the values#
#             of the covariates.                                       #
#                                                                      #
# kernel    : a character. One of "TD" = time dependent,               #
#             "TI" = time-invariant or "LV" = last value.              #
#             This specifies how the kernel is used                    #
#                                                                      #
# lType     : a character. One of "identity", "log", "logistic".       #
#             Specifies the form of the link function.                 #
#                                                                      #
# tt        : If provided, a vector of times at which to evaluate the  #
#             kernel                                                   #
#                                                                      #
#----------------------------------------------------------------------#
#                                                                      #
# Returns scalar value 2.0 * dU %*% U                                  #
#                                                                      #
#----------------------------------------------------------------------#
doptFunc <- function(pars,
                     data.y, 
                     data.x,
                     kernel, 
                     lType,
                     xIs,
                     yIs,
                     nPatients) {

  res <- cmp_duFunc(pars = pars,
                    data.y = data.y, 
                    data.x = data.x,
                    kernel = kernel, 
                    lType = lType,
                    xIs = xIs,
                    yIs = yIs,
                    nPatients = nPatients)

  return(2.0*(res$du %*% res$u))


}
