quantBMAnormal <-
function(alpha, WEIGHTS, MEAN, SD)
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
  lower <- min(MEAN-6*SD)
  upper <- max(MEAN+6*SD)

  if (cdfBMAnormal(lower, WEIGHTS, MEAN, SD, 0) > alpha) return(NA)
  if (cdfBMAnormal(upper, WEIGHTS, MEAN, SD, 0) < alpha) return(NA)

  z <- uniroot(cdfBMAnormal, lower = lower, upper = upper,
               WEIGHTS=WEIGHTS, MEAN=MEAN, SD=SD, offset = alpha)

  z$root
}

