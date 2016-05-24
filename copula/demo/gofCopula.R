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

##>>> NOTA BENE must contain exactly the \dontrun{} part of
## ../man/gofCopula.Rd
## ===================

## A two-dimensional data example ----------------------------------
x <- rCopula(200, claytonCopula(3))

(tau. <- cor(x, method="kendall")[1,2]) # around 0.5 -- 0.6
## Does the Gumbel family seem to be a good choice?
(thG <- iTau(gumbelCopula(), tau.)) # 3.02
gofCopula(gumbelCopula(thG), x)
# SnC: really s..l..o..w.. --- SnB is *EVEN* slower
gofCopula(gumbelCopula(thG), x, method = "SnC")
## What about the Clayton family?
(thC <- iTau(claytonCopula(), tau.)) # 4.05
gofCopula(claytonCopula(thC), x)
gofCopula(claytonCopula(thC), x, method = "AnChisq")

## The same with a different estimation method
gofCopula(gumbelCopula (thG), x, estim.method="itau")
gofCopula(claytonCopula(thC), x, estim.method="itau")


## A three-dimensional example  ------------------------------------
x <- rCopula(200, tCopula(c(0.5, 0.6, 0.7), dim = 3, dispstr = "un"))

## Does the Clayton family seem to be a good choice?
## here starting with the "same" as indepCopula(3) :
(gCi3 <- gumbelCopula(1, dim = 3, use.indepC="FALSE"))
gofCopula(gCi3, x)
## What about the t copula?
t.copula <- tCopula(rep(0, 3), dim = 3, dispstr = "un", df.fixed=TRUE)
## this is *VERY* slow currently %% FIXME ??
gofCopula(t.copula, x)

## The same with a different estimation method
gofCopula(gCi3,     x, estim.method="itau")
gofCopula(t.copula, x, estim.method="itau")

## The same using the multiplier approach
gofCopula(gCi3,     x, simulation="mult")
gofCopula(t.copula, x, simulation="mult")
