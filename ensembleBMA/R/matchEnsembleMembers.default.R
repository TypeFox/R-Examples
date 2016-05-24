matchEnsembleMembers.default <-
function( fit, ensembleData) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 # match ensemble members in ensembleData and fit

 fitMems <- modelMembers(fit)
 
 ensMems <- ensembleMembers(ensembleData)

 if (!is.null(fitMems) && !is.null(ensMems) 
     && length(fitMems) > length(ensMems))
   stop("model fit has more ensemble members than ensemble data")

 WARN <- rep(FALSE,3)
 WARN[1] <- is.null(fitMems) && !is.null(ensMems)
 WARN[2] <- !is.null(fitMems) && is.null(ensMems)
 WARN[3] <- is.null(fitMems) && is.null(ensMems)

 if (any(WARN) && length(fitMems) != length(ensMems))
   stop("model fit and ensemble data differ in ensemble size")

 if (any(WARN))
  warning("cannot check correspondence between model fit and ensemble data members")

 M <- match(fitMems, ensMems, nomatch = 0)
 if (any(!M)) stop("ensembleData is missing a member used in fit")

 M
}

