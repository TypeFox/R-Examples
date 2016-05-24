`modelParameters.ensembleBMAnormal` <-
function(fit, dates = NULL, ...) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#

 dateTable <- dimnames(fit$weights)[[2]]

 if (is.null(dates)) dates <- dateTable

 dates <- sort(unique(as.character(dates)))

 if (length(dates) > length(dateTable)) 
   stop("parameters not available for some dates")

 I <- match( dates, dateTable, nomatch=0)

 if (any(!I) || !length(I)) 
   stop("parameters not available for some dates")

 list(weights = fit$weights[,I], 
      biasCoefs = fit$biasCoefs[,,I], 
      sd = if (is.null(dim(sd))) fit$sd[I] else fit$sd[,I],
      model = "normal")
}

