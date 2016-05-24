library(pcalg)

cat("doExtras:", (doExtras <- pcalg:::doExtras()), "\n")
showProc.time <- local({
    pct <- proc.time()
    function() { ## CPU elapsed __since last called__
	ot <- pct ; pct <<- proc.time()
	cat('Time elapsed: ', (pct - ot)[1:3],'\n')
    }
})
options(warn = 1)# show warnings immediately

load("discr100k.rda")# in directory tests/ i.e., typically *not* installed
## Note: --> ../man/gmD.Rd  is _same_ data, only just the first 10'000 obs.
str(dat)# 100'000 x 5 data frame of integers
sapply(dat, table)
nlev.dat <- sapply(dat, function(v) nlevels(factor(v)))
suffStat <- list(dm = dat, nlev = nlev.dat, adaptDF = TRUE)
dat. <- data.matrix(dat)

showProc.time()

## Check for reasonable error message:
Sys.setlocale("LC_MESSAGES", "C")
(er <- tryCatch(gSquareDis(1,3, S=20, dm=dat.),  error=conditionMessage))
##  "x, y, and S must all be in {1,..,p}, p=5"
stopifnot(grepl("must all be in {1,..,p}", er, fixed=TRUE))


## Tests for   |S| = 0  and  |S| = 1

pv.5 <- c(
    ## Collider
    ##        x y  S
    disCItest(1,2,NULL,suffStat),# 1 : Must be !=0
    disCItest(1,2, 3,  suffStat),# 2 : 0
    disCItest(5,2, 4,  suffStat),# 3 : 0
    disCItest(5,2,NULL,suffStat),# 4 : != 0
    ## Non-collider
    disCItest(4,3, 2,  suffStat),# 5 : !=0
    disCItest(4,3,NULL,suffStat),# 6 : 0
    ## Neighbours
    disCItest(4,2,NULL,suffStat),# 7 : 0
    disCItest(1,3,NULL,suffStat),# 8 : 0
    disCItest(1,5,NULL,suffStat))# 9 : !=0

stopifnot(all.equal(pv.5,
                    c(0.5150866,
                      0,
                      0,
                      0.289990516,
                      ## Non-collider
                      0.818905483,
                      0,
                      # Neighbours
                      0,
                      0,
                      0.774158263), tol = 4e-9))

showProc.time()

if(!doExtras) ## smaller data set
    dat. <- dat.[1:4096,]

##  |S| = 2 :
pv2.5 <- c(
    gSquareDis(3,4, S = c(1,2), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(2,4, S = c(1,3), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(2,3, S = c(1,4), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(2,3, S = c(1,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,4, S = c(2,3), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,3, S = c(2,4), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,3, S = c(2,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,2, S = c(3,4), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,2, S = c(3,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,2, S = c(4,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),

    gSquareDis(3,5, S = c(1,2), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(2,5, S = c(1,3), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(2,5, S = c(1,4), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(2,4, S = c(1,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,5, S = c(2,3), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,5, S = c(2,4), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,4, S = c(2,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,5, S = c(3,4), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,4, S = c(3,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,3, S = c(4,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE))

showProc.time()

pv2.5.EXP <- {
    if(doExtras)
        c(0.98726, 0,       0, 0, 0.77097, 0,       0,       0, 0, 0.98974,
          0.95964, 0.99966, 0, 0, 0.99685, 0.99999, 0.42513, 1, 0.11179, 0)
    else
        c(0.73219, 0, 0.0003349, 0, 0.90758, 0, 0, 0.98106, 0.46011, 1,
          0.93245, 0.99998, 8e-07, 0, 0.99797, 0.99933, 0.86036, 0.99988, 0.84859, 0)
}


(eq52 <- all.equal(pv2.5, pv2.5.EXP, tol = 1e-5))
if(!isTRUE(eq52)) stop("gSquareDis(|S| = 2) not ok:", eq52)


##  |S| = 3 :
pv3.5 <- c(
    gSquareDis(1,2, S = c(3,4,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,3, S = c(2,4,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,4, S = c(2,3,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(1,5, S = c(2,3,4), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(2,3, S = c(1,4,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(2,4, S = c(1,3,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(2,5, S = c(1,3,4), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(3,4, S = c(1,2,5), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(3,5, S = c(1,2,4), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE),
    gSquareDis(4,5, S = c(1,2,3), dm = dat., nlev = nlev.dat, adaptDF = TRUE, verbose=TRUE))
showProc.time()

pv3.5.EXP <- {
    if(doExtras)
        c(0,       0, 0.95861, 1, 0,       0, 0,       0.99204, 1, 0)
    else
        c(0.99999, 0, 0.84316, 1, 0.15977, 0, 0.83432, 0.66329, 1, 0)
}

(eq53 <- all.equal(pv3.5, pv3.5.EXP, tol=6e-6))
if(!isTRUE(eq53)) stop("gSquareDis(|S| = 3) not ok:", eq53)

###  Extend data, and do (few) tests with  |S| in {4, 5, 6} !!
###  ---------------------------------------------------------

## Idee:  1) marginal-table fuer die no-incoming (aka 'root') notes
##        2) bivariate table -->  X_j | X_i  for all edges  i -> j
## ==> randomly generate generate the discrete data from these.

###--> Look at  gR  CRAN task view -- package 'gRain'
###--- can easily construct such objects

##' Generate discrete (integer, in 0:k) random variable, given one parent (in a DAG)
##'
##' @title Generate Discrete Random Variable Given one parent
##' @param pa.var
##' @param Pr Markov / conditional probabability matrix Pr[i,j] = Prob{newX = i, pa.X = j}
##' @return a discrete random variable \code{newX} of the same length as \code{pa.var}
##' @author Martin Maechler
rDiscr1 <- function(pa.var, Pr, zero.based = TRUE) {
    stopifnot(pa.var %% 1 == 0, # <- fast test for integer valued
              is.matrix(Pr), is.numeric(Pr),
              (nlev <- ncol(Pr)) >= 1, ## #{levels} of the (discrete) parent variable
              nlev == nlevels(pa.f <- factor(pa.var)))

    pa.n <- as.integer(pa.f)## values in {1, 2,.., nlev}
    X <- integer(n <- length(pa.n))
    FF <-
        if(zero.based)
             function(.) which(. == 1L) - 1L ## {1,2,3} - 1  ==> {0,1,2}
        else function(.) which(. == 1L)
    for(j in seq_len(nlev)) {
        is.j <- pa.n == j
        nj <- sum(is.j)
        X[is.j] <- apply(rmultinom(nj, size = 1, prob = Pr[, j]),
                         2L, FF)
    }
    X
}

nrow(dat.) # 4096 or 100'000  depending on  doExtras

## Currently, already |S| = 4 is so slow, we need
## smaller data set:
dat. <- dat.[ seq_len((if(doExtras) 8192 else 1024)), ]
(n <- nrow(dat.))

set.seed(102)
P6.5 <- cbind(c(1,5,2),
              c(5,2,1))/8 ## P6.[i,j] :=  P[X6 = i | X5 = j]:
X6 <- rDiscr1(dat.[,"X5"], Pr = P6.5)
P7.3 <- cbind(c(1,5,2),
              c(5,2,1),
              c(1,4,3) )/8
X7 <- rDiscr1(dat.[,"X3"], Pr = P7.3)
X8 <- rDiscr1(       X7  , Pr = P7.3)
X9 <- apply(rmultinom(n, size=1, prob = c(1,2,3,4,6)/16),
            2L, function(.) which(. == 1L) - 1L) # independent
P10.9 <- cbind(c(1,7),
               c(6,2),
               c(3,5),
               c(4,4),
               c(2,6))/8
X10 <- rDiscr1(      X9,  Pr = P10.9)
X11 <- rDiscr1(      X9,  Pr = P10.9)
dat10 <- cbind(dat., X6, X7, X8, X9, X10, X11)
(nlev10 <- sapply(as.data.frame(dat10), function(v) nlevels(factor(v))))
showProc.time()

## Testing only (not good for sound inference!):
## setting  n.min = n-1  ==> the tests are also run  for "too small" sample sizes
##          ===========

gSquareDis.. <- function(x,y, S, dm, nlev, verbose = TRUE)
    gSquareDis(x,y, S=S, dm=dm, nlev=nlev,
               n.min=nrow(dm)-1L, adaptDF=TRUE, verbose=verbose)

##  |S| = 4 :
pv4.10 <- c(
    gSquareDis..(1,2, S = 3:6     , dm= dat10, nlev=nlev10),
    gSquareDis..(1,3, S = 4:7     , dm= dat10, nlev=nlev10),
    gSquareDis..(1,4, S = 5:8     , dm= dat10, nlev=nlev10),
    gSquareDis..(1,5, S = 6:9     , dm= dat10, nlev=nlev10),
    gSquareDis..(2,3, S = 4:7     , dm= dat10, nlev=nlev10),
    gSquareDis..(2,4, S = c(3,5:7), dm= dat10, nlev=nlev10),
    gSquareDis..(2,5, S = c(4,6:8), dm= dat10, nlev=nlev10),
    gSquareDis..(2,5, S = c(4,7:9), dm= dat10, nlev=nlev10),
    gSquareDis..(2,5, S = 6:9     , dm= dat10, nlev=nlev10),
    gSquareDis..(3,4, S = 5:8     , dm= dat10, nlev=nlev10),
    gSquareDis..(3,4, S = 6:9     , dm= dat10, nlev=nlev10),
    gSquareDis..(3,5, S = 6:9     , dm= dat10, nlev=nlev10),
    gSquareDis..(4,5, S = 6:9     , dm= dat10, nlev=nlev10),
    gSquareDis..(4,5, S = c(1:3,6), dm= dat10, nlev=nlev10),
    gSquareDis..(4,5, S = c(1:3,7), dm= dat10, nlev=nlev10),
    gSquareDis..(4,5, S = c(1:3,8), dm= dat10, nlev=nlev10))
showProc.time()
dput(signif(zapsmall(pv4.10), 4))

pv4.10.EXP <- {
  if(doExtras)
    c(1,      0, 0.9757, 1, 0.9981, 0,    1,1,1, 0.6949, 0.7201, 1,   0,    0,0,0)
  else
    c(0.9999, 0, 1e-07,  1, 1, 0.0001146, 1,1,1, 5.9e-06,  0,    1, 0.2205, 0,0,0)
}

(eq4.10 <- all.equal(pv4.10, pv4.10.EXP, tol=5e-5))
if(!isTRUE(eq4.10)) stop("gSquareDis(|S| = 4) not ok:", eq4.10)

##  |S| = 5 : --- verbose=2 gives many messages
##  --------      Adding a new combination of parents at sample  <n>"
pv5.10 <- c(
    gSquareDis..(1,2, S = 3:7     , dm= dat10, nlev=nlev10),
    gSquareDis..(1,3, S = 4:8     , dm= dat10, nlev=nlev10),
    gSquareDis..(1,4, S = 5:9     , dm= dat10, nlev=nlev10),
    gSquareDis..(1,5, S = 6:10    , dm= dat10, nlev=nlev10),
    gSquareDis..(2,3, S = 4:8     , dm= dat10, nlev=nlev10),
    gSquareDis..(2,4, S = c(3,5:8), dm= dat10, nlev=nlev10),
    gSquareDis..(2,5, S = c(4,6:9), dm= dat10, nlev=nlev10),
    gSquareDis..(3,4, S = 5:9     , dm= dat10, nlev=nlev10),
    gSquareDis..(3,4, S = 6:10    , dm= dat10, nlev=nlev10),
    gSquareDis..(3,5, S = 5:9     , dm= dat10, nlev=nlev10),
    gSquareDis..(4,5, S = 6:10    , dm= dat10, nlev=nlev10))
showProc.time()
dput(signif(zapsmall(pv5.10), 5))

pv5.10.EXP <- {
  if(doExtras)
    c(1, 0,      0.037821, 1, 1, 0,       1, 0.0005184, 0.0013204, 1, 0)
  else
    c(1, 0.10938, 0.94484, 1, 1, 0.79863, 1, 0.99495,   0.68719,   1, 1)
}

(eq5.10 <- all.equal(pv5.10, pv5.10.EXP, tol= 5e5))
if(!isTRUE(eq5.10)) stop("gSquareDis(|S| = 5) not ok:", eq5.10)

##  |S| = 6 : --- (see |S|=5):

pv6.10 <- c(
    gSquareDis..(1,2, S = 3:8       , dm= dat10, nlev=nlev10),
    gSquareDis..(1,3, S = 4:9       , dm= dat10, nlev=nlev10),
    gSquareDis..(1,4, S = 5:10      , dm= dat10, nlev=nlev10),
    gSquareDis..(1,5, S = c(2:4,6:9), dm= dat10, nlev=nlev10),
    gSquareDis..(2,3, S = 4:9       , dm= dat10, nlev=nlev10),
    gSquareDis..(2,4, S = 5:10      , dm= dat10, nlev=nlev10),
    gSquareDis..(2,5, S = c(4,6:9)  , dm= dat10, nlev=nlev10),
    gSquareDis..(3,4, S = 5:10      , dm= dat10, nlev=nlev10),
    gSquareDis..(3,4, S = c(1:2,5:8), dm= dat10, nlev=nlev10),
    gSquareDis..(3,5, S = c(1:2,6:9), dm= dat10, nlev=nlev10),
    gSquareDis..(4,5, S = c(1:3,7:9), dm= dat10, nlev=nlev10))
showProc.time()
dput(signif(zapsmall(pv6.10), 5))

pv6.10.EXP <- {
  if(doExtras)
    c(1, 0.023908, 0.000208, 1, 1, 0.94747, 1, 0.0002593, 5.1e-06, 1, 0)
  else
    c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
}

(eq6.10 <- all.equal(pv6.10, pv6.10.EXP, tol= 5e5))
if(!isTRUE(eq6.10)) stop("gSquareDis(|S| = 5) not ok:", eq6.10)
