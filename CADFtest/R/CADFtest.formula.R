CADFtest.formula <- function(model, X=NULL, type=c("trend", "drift", "none"), 
                             data=list(), max.lag.y=1, min.lag.X=0, max.lag.X=0, dname=NULL,
                             criterion=c("none", "BIC", "AIC", "HQC", "MAIC"), ...)
{
# Author:       Claudio Lupi
# This version: December 12, 2008
# This function is an interface to function CADFtest.default that computes Hansen's (1995) CADF test.
# Reference:
# @ARTICLE{,
#   author = {Hansen, Bruce E.},
#   title = {Rethinking the Univariate Approach to Unit Root Testing: {U}sing Covariates to Increase Power},
#   journal = {Econometric Theory},
#   year = {1995},
#   volume = {11},
#   pages = {1148--1171},
#   number = {5},
# }
#
# Arguments:
# model:     a formula. If model is a vector or a time series then the 
#            standard ADF test is performed on the series described by model. If a CADF test is desired, then
#            model should specified in the form z0 ~ z1 + z2 + ... where z0 is the variable to be tested, 
#            while z1 and z2 are the stationary covariates to be used in the test. Note that the model is stylized
#            and all the variables are in levels. It is NOT the model equation on which the test is based.
#            However, the covariates must be STATIONARY. 
# type:      it specifies if the underlying model must be with linear trend ("trend", the default), 
#			 with constant ("drift") or without constant ("none").
# max.lag.y: it specifies the number of lags of the dependent (\Delta y_t).
# min.lag.X: it specifies the maximum lead of the covariates (it must be negative or zero).
# max.lag.X: it specifies the maximum lag of the covariates (it must be positive or zero).
#            then the test is performed using the model that minimizes the selection citerion defined in
#            'criterion'. In this case, the max e min orders serve as upper and lower bounds in the model
#            selection.
# criterion: it can be either "none", "BIC", "AIC" or "HQC". If criterion="none", no automatic model selection
#            is performed. Otherwise, automatic model selection is performed using the specified 
#            criterion.           
#
# The procedure to compute the CADF test p-value is proposed in Costantini et al. (2007). Please cite the paper
# when you use the present function.
# Reference:
# @TECHREPORT{,
#   author = {Costantini, Mauro and Lupi, Claudio and Popp, Stephan},
#   title = {A Panel-{CADF} Test for Unit Roots},
#   institution = {University of Molise},
#   year = {2007},
#   type = {Economics \& Statistics Discussion Paper},
#   number = {39/07},
#   url = {http://econpapers.repec.org/paper/molecsdps/esdp07039.htm}
# }

if (is.null(dname)){dname <- deparse(substitute(model))}

if ((model[3]==".()")|(model[3]=="1()"))
{
  model[3] <- 1
  mf    <- model.frame(model, data=data)
  y     <- model.response(mf)
  X     <- NULL
}
else
{
  mf    <- model.frame(model, data=data)
  y     <- model.response(mf)
  X     <- model.matrix(model, data=data)
  X     <- X[,2:dim(X)[2]]
}

call <- match.call(CADFtest)

test.results <- CADFtest.default(model=y, X=X, type=type, max.lag.y=max.lag.y, 
                                 data=data, min.lag.X=min.lag.X, max.lag.X=max.lag.X, dname=dname,
                                 criterion=criterion, ...)

test.results$call <- call

return(test.results)
}

