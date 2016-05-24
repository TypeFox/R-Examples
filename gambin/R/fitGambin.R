fitGambin <-
function(abundances, subsample = 0)
{
  # for the GamBin function, all the abundances are binned into octaves (on log base 2)

  # test for NA's in the data (and that it is a numeric vector)
  Dataname <- deparse(substitute(abundances))
  
  if(is.vector(abundances) & is.numeric(abundances))
  {
  	mydata <- create_octaves(abundances, subsample)
  } else {
  	if(is.data.frame(abundances))
  	{
  		names(abundances) <- tolower(names(abundances))
  		if(!("species" %in% names(abundances) & "octave" %in% names(abundances))) stop("abundances must be a numeric vector or a data.frame created by create_octaves")	
  		mydata <- abundances[c("octave", "species")]
  	} else stop("abundances must be a numeric vector or a data.frame created by create_octaves")	
  }
  
  val <- optimise(.logLik_gamBin, interval = c(0,30), mydata = mydata)
  res <- list()
  
  res$Alpha <- val[[1]]
  if(res$Alpha == 30) 
  {
    warning("Alpha could not be estimated")
    res$Alpha <- Inf
  }
  
  logLik <- -val[[2]]
  attr(logLik, "df") <- 1
  attr(logLik, "nobs") <- nrow(mydata)
  class(logLik) <- "logLik"
  res$logLik <- logLik
  
  res$fitted.values <- gambin_exp(res$Alpha, max(mydata$octave), sum(mydata$species))


  res$Data <- mydata
  if (subsample == 0)
    res$Dataname <- Dataname else
    res$Dataname <- paste(subsample, "individuals sampled from", Dataname)
  if(length(res$Dataname) > 1) res$Dataname <- "unnamed abundance vector"
  
  res$MaxOctave <- max(mydata$octave)
  
  res$coefficients <- c(Alpha = res$Alpha, MaxOctave = res$MaxOctave)
  
  attr(res, "nobs") <- nrow(mydata)
  class(res) <- "gambin"
  res
}
