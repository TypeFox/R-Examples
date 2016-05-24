################################################################################
# manylm.influence: a function that returns the diagonal(s) of the hat matrix  #
# (if the model is a manyglm, hat is a matrix and each column is the           #
# diag(hat matrix) * dispersion of the respective column of y,                 #
# the change in the estimated coefficients which results when the i-th case    #
# is dropped from the regression,                                              #
# a vector whose i-th element contains the estimate of the residual standard   #
# deviation obtained when the i-th case is dropped from the regression.        #
# and a vector of weighted (or for class glm rather deviance) residuals.       #
#                                                                              #
# Modified by Alice:                                                           #
#   Stop using .Fortran("lminfl"...) in non-base packages as required by CRAN  #
# Instead, estimate multiple attributes manually                               # 
################################################################################

manylm.influence <- function (model, do.coef = TRUE) {

    res$wt.res <- as.matrix(weighted.residuals(model))
    res$hat <- diag(model$hat.X)
    warning("Not implemented for manylm objects")

    return(res)
}


