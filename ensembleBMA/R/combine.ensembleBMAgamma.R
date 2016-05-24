combine.ensembleBMAgamma <-
function (x, y, ...) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#

 if (length(list(...)) == 0) {
   z <- list( x, y)
 }
 else {
   z <- c(list( x, y), list(...))
 }

 if (!all(sapply(z, function(x) class(x)[1]) == "ensembleBMAgamma")) 
      stop("merge components differ in class")

 if (length(power <- unique(sapply(z, function(x) x$power))) != 1) 
      stop("merge components differ in data tranformation")

 core <- combineCORE(z)
 DATES <- names(core$nIter) 
 nam <- core$sortedMemberNames

 n <- length(core$notdup)
 prob0coefs <- array( NA, c( dim(x$prob0coefs)[1:2], n))
 biasCoefs <- array( NA, c( dim(x$biasCoefs)[1:2], n))
 varCoefs <- matrix( NA, nrow(x$varCoefs), n)

 j <- 0
 for (i in 1:length(z)) {
    l <- dim(z[[i]]$weights)[2]
    prob0coefs[,,j+(1:l)] <- z[[i]]$prob0coefs[,nam,]
    biasCoefs[,,j+(1:l)] <- z[[i]]$biasCoefs[,nam,]
    varCoefs[,j+(1:l)] <- z[[i]]$varCoefs
    j <- j + l 
 }

 prob0coefs <- prob0coefs[,,core$notdup]
 biasCoefs <- biasCoefs[,,core$notdup]
 varCoefs <- varCoefs[,core$notdup]

 dimnames(prob0coefs) <- list( NULL, nam, DATES)
 dimnames(biasCoefs) <- list( NULL, nam, DATES)
 dimnames(varCoefs) <- list( NULL, DATES)

# use member ordering as in x 
 nam <- dimnames(x$weights)[[1]]
 prob0coefs <- prob0coefs[,nam,]
 biasCoefs <- biasCoefs[,nam,]

 structure(list( training = list(days = days, lag = lag, table = table), 
          prob0coefs = prob0coefs, biasCoefs = biasCoefs, varCoefs = varCoefs,
          weights = core$weights, nIter = core$nIter, 
          exchangeable = core$exch, power = power), 
                 forecastHour = core$forecastHour,
                 initializationTime = core$initializationTime, 
                 call = lapply( z, function(x) attr(x,"call")),
                 class = class(x))
}

