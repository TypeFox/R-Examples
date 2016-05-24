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

source(system.file("Rsource", "utils.R", package="copula", mustWork=TRUE))

set.seed(101)
X <- rstable1(1e4, alpha=.9999, beta=1, gamma= .25, delta=1)
summary(X)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##   1592    1592    1593    1594    1593    2092
stopifnot(all(X > 1500))
X <- rstable1(1e4, alpha=.999999, beta=1, gamma= .25, delta=1)
stopifnot(150000 < X, X < 170000)
showProc.time()

r1R <- copula:::rstable1R
r1C <- copula:::rstable1C
set.seed(123)

## speed is *very* similar (!):
showSys.time(Z  <- r1R(10000, .2, beta=1, gamma= 45))
showSys.time(Z. <- r1C(10000, .2, beta=1, gamma= 45))

ks.test(Z, Z.) # p-value ~ 0.50  they "are the same"

if(!dev.interactive(orNone=TRUE)) pdf("rstable-ex.pdf")
qqplot(Z, Z., log = "xy") # looks nice
acf(log(Z))  # "nice"
acf(log(Z.)) # ditto

showProc.time()

### Ensure basic properties of  rCopula() ---------------------

pkg <- "package:copula"
## all "acopula"s :
aCops <- sapply(ls(pkg, pattern="^cop[A-Z]"),
                get, pkg, simplify=FALSE)
stopifnot(sapply(aCops, is, "acopula"))

taus <- 1/ c(16, 8, 5, 4) # < 1/3 which is max{ AMH }

aCt2 <- unlist(lapply(taus, function(TAU) {
    lapply(aCops, function(COP) {
        th <- iTau(COP, TAU); onacopulaL(COP, list(th, 1:2))})}))

## Now the IJ-copulas:
(csubc <- names(getClass("copula")@subclasses))
## and only keep the "actual" (i.e. non-virtual) ones:
(csubc <- Filter(Negate(isVirtualClass), csubc))
## schlatherCopula() e.g. has no constructor:
(cfnms <- intersect(csubc, ls(pkg, pattern = "[A-Za-z]+Copula$")))
stopList <- c("asymCopula",  "asymExplicitCopula", "indepCopula")
cfnms <- cfnms[is.na(match(cfnms, stopList))]
str(cfn <- sapply(cfnms, get, pkg, simplify=FALSE))
str(th.25 <- lapply(cfn, function(F) iTau(F(), 0.25)))
Ct2 <- unlist(lapply(taus, function(TAU) {
    lapply(cfn, function(cFn) { th <- iTau(cFn(), TAU); cFn(th)})}))

## A list of "fully specified" (parameter, dimension) copula models :
length(cops <- c(aCt2, Ct2))
## 72 of them ...

set.seed(17)
Uc <- lapply(cops, rCopula, n = 1024)

## 1) Check the tau's

## TODO


## 2) Check the margins: all must be uniform
ad.test <- ADGofTest::ad.test
pp <- sapply(Uc, function(U2)
             apply(U2, 2,
                   function(u) ad.test(u)$p.value))

## The P-values should simply be uniform in [0,1]:
hh <- hist(c(pp), breaks = (0:20)/20)## should "look" uniform
## and hence their test should typically *not* be significant
(ad.pp <- ad.test(c(pp)))
stopifnot(hh$counts[1] < 15, # for the above seed, it is 13 , "somewhat high"
          ad.pp$p.value > 0.10)
