cdf.fitMOSnormal <-
function(fit, ensembleData, values, dates = NULL, ...)
{

 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

## remove instances missing all forecasts or dates

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
## M <- M | is.na(ensembleDates(ensembleData))
 ensembleData <- ensembleData[!M,]

 nObs <- nrow(ensembleData)

 if (!is.null(dates)) warning("dates ignored")

 CDF <- matrix(NA, nObs, length(values))
 dimnames(CDF) <- list(ensembleObsLabels(ensembleData),as.character(values))
 Mu <- rep(NA,nObs)
 Sig <- rep(NA,nObs)

 ensembleData <- ensembleForecasts(ensembleData)
 x <- c(fit$a,fit$B)
 A <- cbind(rep(1,nObs),ensembleData)

 S.sq <- apply(ensembleData,1,var)
 Mu <- A%*%x
 Sig <- sqrt(rep(fit$c,nObs) + rep(fit$d,nObs)*S.sq)

 CDF <- pnorm(values, mean = Mu, sd = Sig)
 CDF
}

