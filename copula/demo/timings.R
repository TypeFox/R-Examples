## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


require(copula)

### Results of using nacFrail.time() ############################################

set.seed(1) # set seed

n <- 10000
taus <- c(0.05,(1:9)/10,0.95)


### AMH

nacFrail.time(n,"AMH",taus[taus < 0.33])

##             inner tau
## outer tau   0.10 0.20 0.30
##      0.05 1    4    4    4
##      0.10 2   NA    4    4
##      0.20 2   NA   NA    5
##      0.30 3   NA   NA   NA

## conclusion:
## V0:  uniformly fast
## V01: uniformly fast


### Clayton

nacFrail.time(n,"Clayton",taus)

##             inner tau
## outer tau   0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95
##      0.05 2   38   44   55   62   56   58   52   47   38   34
##      0.10 2   NA   53   58   57   58   63   61   58   55   51
##      0.20 2   NA   NA   42   41   37   43   44   44   44   42
##      0.30 3   NA   NA   NA   28   28   29   29   30   26   30
##      0.40 3   NA   NA   NA   NA   21   23   22   19   22   23
##      0.50 2   NA   NA   NA   NA   NA   17   17   18   18   19
##      0.60 3   NA   NA   NA   NA   NA   NA   15   15   16   15
##      0.70 2   NA   NA   NA   NA   NA   NA   NA   13   14   13
##      0.80 3   NA   NA   NA   NA   NA   NA   NA   NA   12   12
##      0.90 3   NA   NA   NA   NA   NA   NA   NA   NA   NA   11
##      0.95 2   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA

## conclusion:
## V0:  uniformly fast
## V01: faster for larger parameter


### Frank

nacFrail.time(n,"Frank",taus)

##             inner tau
## outer tau   0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95
##      0.05 1   12   14   13   13   13   13   13   13   12   11
##      0.10 2   NA   21   23   23   24   24   23   23   17   16
##      0.20 3   NA   NA   39   53   62   66   71   74   21   20
##      0.30 3   NA   NA   NA   91  133  150  184  201   40   39
##      0.40 2   NA   NA   NA   NA  175  296  369  389   92   92
##      0.50 3   NA   NA   NA   NA   NA  759 1471 2002  310  318
##      0.60 2   NA   NA   NA   NA   NA   NA 3493 5736 1864 1996
##      0.70 3   NA   NA   NA   NA   NA   NA   NA 2646 4266 4904
##      0.80 5   NA   NA   NA   NA   NA   NA   NA   NA 2126 2900
##      0.90 3   NA   NA   NA   NA   NA   NA   NA   NA   NA  986
##      0.95 3   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA

## conclusion:
## V0:  uniformly fast
## V01: depending on theta0, theta1, and approx (due to rej)


### Gumbel

nacFrail.time(n,"Gumbel",taus)

##             inner tau
## outer tau   0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95
##      0.05 8    9    9    9    9    9    9    9    9    9    9
##      0.10 8   NA    9    9    9    9    9    9    9    9   10
##      0.20 8   NA   NA    9    9    8    5   10    9    9    9
##      0.30 7   NA   NA   NA    9    8   10   10    9   10    9
##      0.40 8   NA   NA   NA   NA    9    9    9    8    7    9
##      0.50 5   NA   NA   NA   NA   NA    9    9    9   10   10
##      0.60 8   NA   NA   NA   NA   NA   NA    9    9    9   10
##      0.70 7   NA   NA   NA   NA   NA   NA   NA    9   10   10
##      0.80 8   NA   NA   NA   NA   NA   NA   NA   NA    6    9
##      0.90 7   NA   NA   NA   NA   NA   NA   NA   NA   NA    8
##      0.95 8   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA

## conclusion:
## V0:  uniformly fast
## V01: uniformly fast


### Joe

nacFrail.time(n,"Joe",taus)

##             inner tau
## outer tau   0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95
##      0.05 1    2    3    5    6    7    8    9    9   10    9
##      0.10 2   NA    7   12   15   19   21   24   25   26   25
##      0.20 3   NA   NA   22   36   49   60   69   75   78   77
##      0.30 4   NA   NA   NA   53   89  121  147  167  178  177
##      0.40 5   NA   NA   NA   NA  154  264  359  433  481  483
##      0.50 4   NA   NA   NA   NA   NA  330  570  768  910  934
##      0.60 6   NA   NA   NA   NA   NA   NA  650 1129 1493 1598
##      0.70 5   NA   NA   NA   NA   NA   NA   NA 1031 1774 2034
##      0.80 6   NA   NA   NA   NA   NA   NA   NA   NA 1680 2280
##      0.90 6   NA   NA   NA   NA   NA   NA   NA   NA   NA 1298
##      0.95 6   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA

## conclusion:
## V0:  (almost) uniformly fast
## V01: increasing in theta0 and theta1
## V01 depends on theta0 and theta1 [but is bounded due to approx-parameter]
