`quantBMAgamma0` <-
function(alpha, WEIGHTS, MEAN, VAR, PROB0)
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#

 # if the probability of zero is greater than the desired alpha
 # then the quantile is zero

  if (sum(WEIGHTS*PROB0) > alpha) return(0)

# Initialize: Find lower and upper bounds
 
  lower <- 0
  upper <- max(MEAN+6*sqrt(VAR))

  if (cdfBMAgamma0(lower, WEIGHTS, MEAN, VAR, PROB0, 0) > alpha) return(NA)
  if (cdfBMAgamma0(upper, WEIGHTS, MEAN, VAR, PROB0, 0) < alpha) return(NA)

  z <- uniroot(cdfBMAgamma0, lower = lower, upper = upper,
           WEIGHTS=WEIGHTS, MEAN=MEAN, VAR=VAR, PROB0 = PROB0, offset = alpha)

# print(c(alpha, z$root,abs(z$f.root)))

  z$root
}

