`cdfBMAnormal` <-
function (x, WEIGHTS, MEAN, SD, offset = 0)
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
  sum(WEIGHTS*pnorm(x, mean = MEAN, sd = SD)) - offset
}

