matchEnsembleMembers.ensembleMOSnormal <-
function(fit, ensembleData)
{
 # match ensemble members in ensembleData and fit

 if (!is.null(dim(fit$B))) {
   fitMems <- dimnames(fit$B)[[1]]
 }
 else {
   fitMems <- names(fit$B)
 }

 ensMems <- ensembleMemberLabels(ensembleData)

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

