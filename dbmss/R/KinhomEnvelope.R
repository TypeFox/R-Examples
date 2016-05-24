KinhomEnvelope <-
function(X, r = NULL, NumberOfSimulations = 100, Alpha = 0.05, 
         ReferenceType = "", lambda = NULL, SimulationType = "RandomPosition", Global = FALSE) {

  CheckdbmssArguments()
  
  # Estimate intensity if it has not been provided.
  if (is.null(lambda)) {
    if (ReferenceType == "") {
      X.reduced <- X
    } else {
      X.reduced <- X[X$marks$PointType==ReferenceType]
    }
    lambda <- density.ppp(X.reduced, sigma=bw.diggle(X.reduced))
  }
  
  # Choose the null hypothesis
  SimulatedPP <- switch (SimulationType,
                         RandomPosition = expression(rpoispp(lambda)),
                         RandomLocation = expression(rRandomLocation(X, CheckArguments = FALSE)),
                         RandomLabeling = expression(rRandomLabeling(X, CheckArguments = FALSE)),
                         PopulationIndependence = expression(rPopulationIndependenceM(X, ReferenceType = ReferenceType, CheckArguments = FALSE))
                         )
  if (is.null(SimulatedPP))
    stop(paste("The null hypothesis", sQuote(SimulationType), "has not been recognized."))
  # local envelope, keep extreme values for lo and hi (nrank=1)
  Envelope <- envelope(X, fun=Kinhomhat, nsim=NumberOfSimulations, nrank=1,
                       r=r, ReferenceType=ReferenceType, lambda=lambda, 
                       CheckArguments = FALSE,
                       simulate=SimulatedPP, savefuns=TRUE
                       )
  attr(Envelope, "einfo")$H0 <- switch (SimulationType,
                                        RandomPosition = "Random Position",
                                        RandomLocation = "Random Location",
                                        )
  # Calculate confidence intervals
  Envelope <- FillEnvelope(Envelope, Alpha, Global)
  # Return the envelope
  return (Envelope)
}
