CRPS.default <-
function(fit, ensembleData, dates=NULL, nSamples=NULL, seed=NULL, ...) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 mc <- match.call()   
 mc[[1]] <- as.name("crps")
 colMeans( eval(mc,parent.frame()), na.rm=TRUE)
}

