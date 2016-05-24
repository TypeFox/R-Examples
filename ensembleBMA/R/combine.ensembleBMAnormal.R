combine.ensembleBMAnormal <-
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

 if (!all(sapply(z, function(x) class(x)[1]) == "ensembleBMAnormal")) 
      stop("merge components differ in class")

 core <- combineCORE(z)
 DATES <- names(core$nIter)
 nam <- core$sortedMemberNames

 sd <- unlist(lapply( z, function(x) x$sd))[core$notdup][DATES]
 
 biasCoefs <- array( NA, c( dim(x$biasCoefs)[1:2], length(core$notdup)))

 j <- 0
 for (i in 1:length(z)) {
    l <- dim(z[[i]]$weights)[2]
    biasCoefs[,,j+(1:l)] <- z[[i]]$biasCoefs[,nam,]
    j <- j + l 
 }

 biasCoefs <- biasCoefs[,,core$notdup]

 dimnames(biasCoefs) <- list( NULL, nam, DATES)

# use member ordering as in x 
 nam <- dimnames(x$weights)[[1]]
 biasCoefs <- biasCoefs[,nam,]

 structure(list( training = core$training,
                 biasCoefs = biasCoefs, sd = sd,
       weights = core$weights, nIter = core$nIter, exchangeable = core$exch), 
                 forecastHour = core$forecastHour,
                 initializationTime = core$initializationTime, 
                 call = lapply( z, function(x) attr(x,"call")),
                 class = class(x))
}

