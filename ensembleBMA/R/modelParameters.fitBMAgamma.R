`modelParameters.fitBMAgamma` <-
function(fit, ...) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 list(weights = fit$weights, 
      popCoefs = fit$popCoefs,
      biasCoefs = fit$biasCoefs, 
      varCoefs = fit$varCoefs,
      power = fit$power,
      model = "gamma")
}

