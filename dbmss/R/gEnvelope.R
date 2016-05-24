gEnvelope <-
function(X, r = NULL, NumberOfSimulations = 100, Alpha = 0.05, 
         ReferenceType = "", NeighborType = "", SimulationType = "RandomPosition", Global = FALSE) {

  CheckdbmssArguments()
  
  # Choose the null hypothesis
  SimulatedPP <- switch (SimulationType,
                         RandomPosition = expression(rRandomPositionK(X, CheckArguments = FALSE)),
                         RandomLabeling = expression(rRandomLabeling(X, CheckArguments = FALSE)),
                         PopulationIndependence = expression(rPopulationIndependenceK(X, ReferenceType, NeighborType, CheckArguments = FALSE))
                         )
  if (is.null(SimulatedPP))
    stop(paste("The null hypothesis", sQuote(SimulationType), "has not been recognized."))
  # local envelope, keep extreme values for lo and hi (nrank=1)
  Envelope <- envelope(X, fun=ghat, nsim=NumberOfSimulations, nrank=1,
                       r=r, ReferenceType=ReferenceType, NeighborType=NeighborType, 
                       CheckArguments = FALSE,
                       simulate=SimulatedPP, savefuns=TRUE
                       )
  attr(Envelope, "einfo")$H0 <- switch (SimulationType,
                                        RandomPosition = "Random Position",
                                        RandomLabeling = "Random Labeling",
                                        PopulationIndependence = "Population Independence"
                                        )
  # Calculate confidence intervals
  Envelope <- FillEnvelope(Envelope, Alpha, Global)
  # Return the envelope
  return (Envelope)
}
