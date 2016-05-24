getInitiallbc <-
function(theta, sigmaP)
{
   if (length(theta) != length(sigmaP))
      stop("incompatible dimensions!")
   dat <- data.frame(theta, sigmaP)
   ini <- getInitial(sigmaP ~ SSlbc(theta, b0, b1),
      data = dat)
   return(ini)
}
