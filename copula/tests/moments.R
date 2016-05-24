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

if(getRversion() < "2.15")
paste0 <- function(...) paste(..., sep="")

source(system.file("Rsource", "utils.R", package="copula", mustWork=TRUE))
##-> setPar(), showProc.time() etc
## All non-virtual copula classes:
source(system.file("Rsource", "cops.R", package="copula", mustWork=TRUE))
## --> copcl, copcl., copObs, copBnds,  excl.2 , copO.2, copBnd.2

showProc.time()

copcl ## the classes (incl. 'indepCopula')
str(copObs, max.level=1)# copula objects (w/o 'indepCopula')
## and their parameter bounds:
t(copBnds)

###-------- tau & and inverse ---------------------------------------------------

## currently fails: --- FIXME?: should AMH also warn like the others?
tau.s <- c(-.999, -.1, 0, (1:3)/10, .5, .999)
### give different warnings , but "work" :  { .5 , even 1/3, gives error for AMH FIXME}

##
tau(tevCopula(0)) # 0.05804811
## restricted tau-range works better
tau.s <- c(       -.1, 0, 0.05805, (1:2)/9, 0.3)

names(tau.s) <- paste0("tau=", sub("0[.]", ".", formatC(tau.s)))
tTau <- sapply(tau.s, function(tau)
               sapply(copObs, iTau, tau = tau))
tTau
tTau["joeCopula", "tau=-.1"] <- 1 # tauJoe() "works outside admissible range"

stopifnot(rep(copBnds["min",],ncol(tTau)) <= tTau + 1e-7,
          tTau <= rep(copBnds["max",],ncol(tTau)),
          ## theta and tau are comonotone :
          apply(tTau, 1, diff) >= -1e9)

showProc.time()

tautau <- t(sapply(names(copObs), function(cNam)
		   sapply(tTau[cNam,],
			  function(th) tau(setPar(copObs[[cNam]], th)))))
tautau
showProc.time()

xctTau <- matrix(tau.s, nrow = nrow(tautau), ncol=length(tau.s),
                 byrow=TRUE)
## The absolute errors
errTau <- tautau-xctTau
round(10000*errTau)
## has two NaN .. ok, for now:
errTau["tawnCopula", 1:2] <- 0
## These families do not support tau < 0
errTau[c("gumbelCopula", "joeCopula",
         "galambosCopula", "huslerReissCopula","tevCopula"),
       "tau=-.1"] <- 0
## the tevCopula cannot get a tau = 0 (for now) __FIXME?__
errTau["tevCopula", 2] <- 0
## "fgmCopula" has tau in [-2/9, 2/9] :
errTau["fgmCopula", "tau=.3"] <- 0
stopifnot(max(abs(errTau)) <= 0.00052)# ok for IJ-taus

showProc.time()


###-------- rho & and inverse ---------------------------------------------------

## NB:
##  iRho() method for class "amhCopula" not yet implemented

### give different warnings , but "work" [not using AMH and Joe]:
rho.s <- c(-.999, -.1, 0, (1:3)/9, .5, .9, .999)
names(rho.s) <- paste0("rho=", sub("0[.]", ".", formatC(rho.s)))
tRho <- sapply(rho.s, function(rho)
               sapply(copO.2, iRho, rho = rho))
tRho
warnings()## 16 warnings [2014-05]
## and from now on, show them as they happen:
options(warn = 1)

##--> oops!  clayton [rho=0] is NA __FIXME__
tRho["claytonCopula", "rho=0"] <- 0
## and it has NA also for .999  {but that maybe consider ok}:
tRho["claytonCopula", "rho=.999"] <- 10^100

stopifnot(rep(copBnd.2["min",],ncol(tRho)) <= tRho,
          tRho <= rep(copBnd.2["max",],ncol(tRho)),
          ## theta and rho are comonotone :
          apply(tRho, 1, diff) >= 0)

rhorho <- t(sapply(names(copO.2), function(cNam)
                   sapply(tRho[cNam,],
                          function(th) rho(setPar(copO.2[[cNam]], th)))))
rhorho
showProc.time()

xctRho <- matrix(rho.s, nrow = nrow(rhorho), ncol=length(rho.s),
                 byrow=TRUE)
## The absolute errors
errRho <- rhorho-xctRho
round(10000*errRho)
## the tevCopula cannot get a rho <= 0 (for now) __FIXME?__
errRho["tevCopula", 1:3] <- 0
## These three families do not support rho < 0  (currently):
errRho[c("gumbelCopula", "galambosCopula", "huslerReissCopula", "tawnCopula"),
       c("rho=-.1","rho=-.999")] <- 0
errRho["tawnCopula", rho.s >= 0.9] <- 0
## "fgmCopula" has rho in [-1/3, 1/3] :
errRho["fgmCopula", abs(rho.s) > 1/3] <- 0

stopifnot(max(abs(errRho)) <= 0.00369,
          max(abs(errRho[,rho.s <= 0.9])) <= 0.0002)# ok for IJ-rhos

showProc.time()


