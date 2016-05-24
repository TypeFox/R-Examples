mpsegrat <- function(object)
{
  if (!inherits(object, "mpcross")) stop("Object must be of class mpcross")

  finals <- object$finals
  founders <- object$founders
  nmrk <- ncol(founders)
  chisq <- vector(length=nmrk)
  pval <- vector(length=nmrk)
  badmrk <- vector()

  for (i in 1:nmrk)
  {
    obs <- table(finals[,i])
    exp <- table(founders[,i])/nrow(founders)*(sum(obs))
    chisq[i] <- NA
    pval[i] <- NA

    if (length(exp)<length(obs))
	badmrk <- c(badmrk, i) else {
	if (length(exp)>length(obs)) {
        obs2 <- vector(length=length(exp))
        obs2[match(names(obs), names(exp))] <- obs
        } else obs2 <- obs

        chisq[i] <- sum((obs2-exp)^2/exp)
        pval[i] <- 1-pchisq(chisq[i], length(obs)-1)
	}
   }
   if (length(badmrk)>0)
	cat("Markers ", badmrk, " had values appear in finals which are not in founders.\n They probably have genotyping errors.\n") 

   df <- data.frame(MarkerName=colnames(founders), chisq, pval)
   df
}

