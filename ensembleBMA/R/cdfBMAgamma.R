`cdfBMAgamma` <-
function (x, WEIGHTS, MEAN, VAR, offset = 0)
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 RATE <- MEAN/VAR
 sum(WEIGHTS*pgamma(x,shape=RATE*MEAN,rate=RATE))-offset
}

