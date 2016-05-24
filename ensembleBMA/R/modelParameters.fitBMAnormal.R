`modelParameters.fitBMAnormal` <-
function(fit, ...) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 list(weights = fit$weights, 
      biasCoefs = fit$biasCoefs, 
      sd = fit$sd,
      model = "normal")
}

