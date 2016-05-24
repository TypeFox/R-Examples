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

user <- Sys.getenv("USER")
switch(user,
       "mhofert"={ # account Marius Hofert (local)
           rm(list=objects()) # remove previously loaded objects
           setwd("~/R/MMMH/parallel/") ## ?????
       },
       "maechler"={
           setwd("~/R/MM/Pkg-ex/copula")
           Mlibrary(copula)
       },
       stop("unsupported user"))

options(warn=1) # print warnings as they occur

## load packages
stopifnot(require(copula))

stime <- function() paste(vapply(unclass(as.POSIXlt(Sys.time()))[3:1],
                                 round, 0.12, digits=2), collapse=":")
stime()
## TODO: improve, the above gives e.g.
## [1] "9:19:31.1"

### Simulation  a la  Genest et al (2009) -- of only one case
### -------------------------------------
gofC <- function(copula, H0copula, n, N, method, seed=NULL,
                 verbose=interactive(), heartbeat=0, fileHB, ...) {
    if(!is.null(seed)) set.seed(seed)
    if(heartbeat) { ii <- if(!is.null(seed)) seed else heartbeat
                    writeLines(sprintf("%s:%5d", stime(), ii), con = fileHB)
                }
    T <- system.time(
        pV <- gofCopula(copula, x = rCopula(n, H0copula), N=N,
                        method=method, verbose=verbose, ...)$p.value)
    c(pvalue = pV, T)# or rather list?
}

nmCop <- function(copula, n = 4)
    substr(sub("Copula$", "", as.vector(class(copula))), 1,n)


## As Genest et al.
n <- 150
dim <- 2
N <- 1000
tau <- 0.5

## 1)
method <- "SnC"
H0copF <- claytonCopula
  copF <- normalCopula
## 2)
method <- "SnB"
H0copF <- tCopula
  copF <- gumbelCopula

H0copF <- claytonCopula
  copF <- normalCopula

if(interactive()) {## Testing:
    N <- 16
    n <- 12
    nRep <- 3
} else {
    nRep <- 2000
}##         ===== {MM: guess: ~ 6 hours, on 24 cores, ada-7 }
## Memory growth : VIRT =  400M (RES =  220M) after  1h30
##                 VIRT =  594M (RES =  313M) after  2h30
##                 VIRT =  662M (RES =  481M) after  4h30
##                 VIRT =  713M (RES =  532M) after  5h07
##                 VIRT =  741M (RES =  560M) after  5h37
##                 VIRT =  772M (RES =  591M) after  6h05
##                 VIRT =  832M (RES =  651M) after  6h32
##                 VIRT =  868M (RES =  687M) after  7h03
##                 VIRT =  942M (RES =  761M) after  8h00
##                 VIRT =  981M (RES =  800M) after  8h35
##                 VIRT = 1024M (RES =  842M) after  9h11
## --> Jobs ends in 10h 30 min {but with a batch error, because I edited this file}


##  nRep <- 10000 ---> on 24 cores (ada-7, Aug.2012) starts accumulating memory
##                ---> about 3 GB per node ..-> machine starts swapping!
##                ---> not finished after > 40 hours ==> 'pkill -f exec/R'

## we set seed *INSIDE* set.seed(17)

(thet.0 <- iTau(H0copF(), tau)) # 2
 cop <- copF(dim=dim)
Hcop <- H0copF(thet.0, dim=dim)

stopifnot(require("parallel"))

filHead <- paste0("gof-sim_", method, "_", nmCop(cop), ":", nmCop(Hcop),
                  "_n=", n,
                  "_d=", dim,
                  "_tau=",formatC(tau),
                  "_nrep=",nRep)
(sFile <- paste0(filHead, ".rda"))
fileHB <- file(paste0(filHead, ".heartbeat"))
open(fileHB, open="wt")
## in case it already exists, start re-writing at beginning:
seek(fileHB, 0); truncate(fileHB)


rr <- lapply(1:nRep, function(i)
	       gofC(cop, H0copula = Hcop, n = n, N = N,
		    seed = i, method = method,
                    heartbeat=TRUE, fileHB=fileHB))

               mc.set.seed=FALSE,
               mc.cores = detectCores())

save(rr, n,dim,N,method, cop, Hcop,
     file = sFile)
close(fileHB)
