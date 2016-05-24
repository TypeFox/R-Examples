## ----echo=F--------------------------------------------------------------
### get knitr just the way we like it

knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE,
  error = FALSE,
  tidy = FALSE,
  cache = FALSE
)

## ------------------------------------------------------------------------
library(stats4)
set.seed(17822)
y <- rpois(n = 40, lambda=10)
# Easy one-dimensional MLE:
nLL <- function(lambda) -sum(stats::dpois(y, lambda, log = TRUE))
fit0 <- mle(nLL, start = list(lambda = 5), nobs = NROW(y))

## ------------------------------------------------------------------------
summary(fit0)

## ------------------------------------------------------------------------
logLik(fit0)

## ------------------------------------------------------------------------
y1 <- y[1:20]
y2 <- y[21:27]
y3 <- y[28:40]

## ------------------------------------------------------------------------
library(gmp)
library(homomorpheR)
Site <- R6::R6Class("Site",
                    private = list(
                        ## name of the site
                        name = NA,
                        ## only master has this, NA for workers
                        privkey = NA,
                        ## local data
                        data = NA,
                        ## The next site in the communication: NA for master
                        nextSite = NA,
                        ## is this the master site?
                        iAmMaster = FALSE,
                        ## intermediate result variable
                        intermediateResult = NA
                    ),
                    public = list(
                        ## Common denominator for approximate real arithmetic
                        den = NA,
                        ## The public key; everyone has this
                        pubkey = NA,
                        initialize = function(name, data, den) {
                            private$name <- name
                            private$data <- data
                            self$den <- den
                        },
                        setPublicKey = function(pubkey) {
                            self$pubkey <- pubkey
                        },
                        setPrivateKey = function(privkey) {
                            private$privkey <- privkey
                        },
                        ## Make me master
                        makeMeMaster = function() {
                            private$iAmMaster <- TRUE
                        },
                        ## add neg log lik and forward to next site
                        addNLLAndForward = function(lambda, enc.offset) {
                            if (private$iAmMaster) {
                                ## We are master, so don't forward
                                ## Just store intermediate result and return
                                private$intermediateResult <- enc.offset
                            } else {
                                ## We are workers, so add and forward
                                ## add negative log likelihood and forward result to next site
                                ## Note that offset is encrypted
                                nllValue <- self$nLL(lambda)
                                result.int <- floor(nllValue)
                                result.frac <- nllValue - result.int
                                result.fracnum <- as.bigq(numerator(as.bigq(result.frac) * self$den))
                                pubkey <- self$pubkey
                                enc.result.int <- pubkey$encrypt(result.int)
                                enc.result.fracnum <- pubkey$encrypt(result.fracnum)
                                result <- list(int = pubkey$add(enc.result.int, enc.offset$int),
                                               frac = pubkey$add(enc.result.fracnum, enc.offset$frac))
                                private$nextSite$addNLLAndForward(lambda, enc.offset = result)
                            }
                            ## Return a TRUE result for now.
                            TRUE
                        },
                        ## Set the next site in the communication graph
                        setNextSite = function(nextSite) {
                            private$nextSite <- nextSite
                        },
                        ## The negative log likelihood
                        nLL = function(lambda) {
                            if (private$iAmMaster) {
                                ## We're master, so need to get result from sites
                                ## 1. Generate a random offset and encrypt it
                                pubkey <- self$pubkey
                                offset <- list(int = random.bigz(nBits = 256),
                                               frac = random.bigz(nBits = 256))
                                enc.offset <- list(int = pubkey$encrypt(offset$int),
                                                   frac = pubkey$encrypt(offset$frac))
                                ## 2. Send off to next site
                                throwaway <- private$nextSite$addNLLAndForward(lambda, enc.offset)
                                ## 3. When the call returns, the result will be in
                                ##    the field intermediateResult, so decrypt that.
                                sum <- private$intermediateResult
                                privkey <- private$privkey
                                intResult <- as.double(privkey$decrypt(sum$int) - offset$int)
                                fracResult <- as.double(as.bigq(privkey$decrypt(sum$frac) - offset$frac) / den)
                                intResult + fracResult
                            } else {
                                ## We're worker, so compute local nLL
                                -sum(stats::dpois(private$data, lambda, log = TRUE))
                            }
                        })
                    )

## ------------------------------------------------------------------------
keys <- PaillierKeyPair$new(1024) ## Generate new public and private key.
den <- gmp::as.bigq(2)^256  #Our denominator for rational approximations

## ------------------------------------------------------------------------
site1 <- Site$new(name = "Site 1", data = y1, den = den)
site2 <- Site$new(name = "Site 2", data = y2, den = den)
site3 <- Site$new(name = "Site 3", data = y3, den = den)

## ------------------------------------------------------------------------
## Master has no data!
master <- Site$new(name = "Master", data = c(), den = den)
master$makeMeMaster()

## ------------------------------------------------------------------------
site1$setPublicKey(keys$pubkey)
site2$setPublicKey(keys$pubkey)
site3$setPublicKey(keys$pubkey)
master$setPublicKey(keys$pubkey)

## ------------------------------------------------------------------------
master$setPrivateKey(keys$getPrivateKey())

## ------------------------------------------------------------------------
master$setNextSite(site1)
site1$setNextSite(site2)
site2$setNextSite(site3)
site3$setNextSite(master)

## ------------------------------------------------------------------------
fit1 <- mle(master$nLL, start = list(lambda = 5))

## ------------------------------------------------------------------------
summary(fit1)

## ------------------------------------------------------------------------
logLik(fit1)

