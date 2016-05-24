#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2015.
##   This file is part of the R package rugarch.
##
##   The R package rugarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rugarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

# A number of loss functions for model evaluation

# squared error
loss.se = function(actual, forecast)
{
	(actual - forecast)^2
}
# absolute error
loss.ae = function(actual, forecast)
{
	abs(actual - forecast)
}
# Trading Return
loss.tr = function(actual, forecast)
{
	-sign(forecast) * actual
}
# Directional Accuracy
loss.da = function(actual, forecast)
{
	-sign(forecast) * sign(actual)
}

# asymmetric error
# == quantile error|given coverage level alpha when using predicted quantiles
# also called 'check loss function', 'tick function' or 'lil-lin loss'
loss.ase = function(actual, forecast, alpha)
{
	alpha - ifelse(actual<forecast, 1, 0)*(actual-forecast)
}
# asymmetric quadratic loss
# when a=0.5 == loss.se
loss.asq = function(actual, forecast, alpha)
{
	alpha - ifelse(actual<forecast, 1, 0)*(actual-forecast)^2
}

# linex function of Varian (1975)
# case alpha>0: exponential for e>0, linear for e<0
# case alpha<0: exponential for e<0, linear for e>0
loss.lin = function(actual, forecast, alpha)
{
	e = actual - forecast
	exp(alpha*e) - alpha*e - 1
}

# double linex Granger (1999)
# exponential for all values of alpha
# alpha=beta == symmetric double linex
loss.2lin = function(actual, forecast, alpha, beta)
{
	e = actual - forecast
	exp(alpha*e) + exp(-beta*e) - (alpha-beta)*e - 2
}
