### integrators ###

# generic method for training integration methods
setGeneric("Regression.Integrator.Fit", function(object, X, y, print.level=1) standardGeneric("Regression.Integrator.Fit"))

# base class for integrator configurations
setClass("Regression.Integrator.Config", slots=c(errfun="function"), contains="VIRTUAL")

# base class for output of integrator training - regression
setClass("Regression.Integrator.FitObj", slots = c(config="Regression.Integrator.Config", est="ANY", pred="numeric"), contains = "VIRTUAL")


### regression operations ###

## select operation ##

# generic method for training integration (i) select methods
setGeneric("Regression.Select.Fit", function(object, X, y, print.level=1) standardGeneric("Regression.Select.Fit"))

# base class for select configurations
setClass("Regression.Select.Config", slots=c(errfun="function"), contains="VIRTUAL")

setClassUnion("RegressionSelectPred", c("NULL","numeric","matrix"))

# base class for output of select training - regression
setClass("Regression.Select.FitObj", slots = c(config="Regression.Select.Config", est="ANY", pred="RegressionSelectPred"), contains = "VIRTUAL")

