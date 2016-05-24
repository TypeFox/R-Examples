## ----noname,echo=FALSE---------------------------------------------------
# opts_chunk$set(comment=NA,background='white')

## ----var-----------------------------------------------------------------
library(mAr)
p <- 3
A <- diag(c(.1, .5, .8))
C <- diag(rep(2, 3))

set.seed(100)
chain <- mAr.sim(w = rep(2,p), A = A, C = C, N = 10000)

## ----means---------------------------------------------------------------
 colMeans(chain)

## ----g-------------------------------------------------------------------
g <- function(x)
{
	return(sum(x^2))
}

## ----mcse----------------------------------------------------------------
library(mcmcse)
mcerror_bm <- mcse.multi(x = chain, method =  "bm", 
	size = "sqroot", g = NULL, level = .95, large = FALSE)
mcerror_bart <- mcse.multi(x = chain, method =  "bartlett", 
	size = "cuberoot", g = NULL, level = .95, large = FALSE)
mcerror_tuk <- mcse.multi(x = chain, method =  "tukey", 
	size = "sqroot", g = NULL, level = .95, large = FALSE)

## ----outputvalue---------------------------------------------------------
mcerror_bm$cov

mcerror_bart$cov

mcerror_tuk$cov

rbind(mcerror_bm$est, mcerror_bart$est, mcerror_tuk$est)

c(mcerror_bm$vol, mcerror_bart$vol, mcerror_tuk$vol)

## ----uni-----------------------------------------------------------------
mcse(x = chain[,1], method = "bm", g = NULL)
mcse.mat(x = chain, method = "bm", g = NULL)

## ----sigma_g-------------------------------------------------------------
g
mcerror_g_bm <- mcse.multi(x = chain, g = g)

mcerror_g_bm$cov

mcerror_g_bm$est

## ----confRegion, out.height = '8cm'--------------------------------------
plot(confRegion(mcerror_bm, which = c(1,2), level = .90), type = 'l', asp = 1)
lines(confRegion(mcerror_bart, which = c(1,2), level = .90), col = "red")

## ----comp_region, out.height = '8cm'-------------------------------------
plot(confRegion(mcerror_bm, which = c(1,2), level = .95), type = 'l', asp = 1)
lines(confRegion(mcerror_bm, which = c(1,2), level = .90), col = "red")

## ----minESS--------------------------------------------------------------
# For mu
minESS(p = 3, alpha = .05, eps = .05)

#For mu_g
minESS(p = 1, alpha = .05, eps = .05)

## ----ess-----------------------------------------------------------------
ess(chain)

## ----multiess------------------------------------------------------------
multiESS(chain)
multiESS(chain, covmat = mcerror_bart$cov)

## ----moresamples---------------------------------------------------------
set.seed(100)
chain <- mAr.sim(w = rep(2,p), A = A, C = C, N = 28000)
multiESS(chain)

## ----estvssamp, out.width = '8cm'----------------------------------------
estvssamp(chain[,1])

## ----qqbig---------------------------------------------------------------
p <- 50
A <- diag(seq(.1, .9, length = p))
C <- diag(rep(2, p))

set.seed(100)
chain <- mAr.sim(w = rep(2,p), A = A, C = C, N = 10000)

## ----qq, out.width = '8cm'-----------------------------------------------
mcerror_bm <- mcse.multi(chain, method = "bm")
qqTest(x = chain, covmat = mcerror_bm$cov)

