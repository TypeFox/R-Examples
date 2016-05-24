`cdfBMAgamma0` <-
function (x, WEIGHTS, MEAN, VAR, PROB0, offset = 0)
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 RATE <- MEAN/VAR
 sum(WEIGHTS*(PROB0+(1-PROB0)*pgamma(x,shape=MEAN*RATE,rate=RATE)))-offset
}

