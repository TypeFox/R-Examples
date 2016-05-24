`quantBMAgamma` <-
function(alpha, WEIGHTS, MEAN, VAR)
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#

# Initialize: Find lower and upper bounds
 
  lower <- 0
  upper <- max(MEAN+6*sqrt(VAR))

  if (cdfBMAgamma(lower, WEIGHTS, MEAN, VAR, 0) > alpha) return(NA)
  if (cdfBMAgamma(upper, WEIGHTS, MEAN, VAR, 0) < alpha) return(NA)

  z <- uniroot(cdfBMAgamma, lower = lower, upper = upper,
               WEIGHTS=WEIGHTS, MEAN=MEAN, VAR=VAR, offset = alpha)

# print(c(alpha, z$root,abs(z$f.root)))

  z$root
}

