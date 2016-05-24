mEnvelope <-
function(X, r = NULL, NumberOfSimulations = 100, Alpha = 0.05, 
         ReferenceType, NeighborType = ReferenceType, CaseControl = FALSE, 
         Original = TRUE, Approximate = ifelse(X$n < 10000, 0, 1), Adjust = 1, 
         MaxRange = "ThirdW", SimulationType = "RandomLocation", Global = FALSE) {
  
  CheckdbmssArguments()
  
  # Choose the null hypothesis
  SimulatedPP <- switch (SimulationType,
                         RandomLocation = expression(rRandomLocation(X, CheckArguments = FALSE)),
                         RandomLabeling = expression(rRandomLabelingM(X, CheckArguments = FALSE)),
                         PopulationIndependence = expression(rPopulationIndependenceM(X, ReferenceType, CheckArguments = FALSE))
                        )
  if (is.null(SimulatedPP))
    stop(paste("The null hypothesis", sQuote(SimulationType), "has not been recognized."))
  # local envelope, keep extreme values for lo and hi (nrank=1)
  Envelope <- envelope(X, fun=mhat, nsim=NumberOfSimulations, nrank=1,
                       r=r, ReferenceType=ReferenceType, NeighborType=NeighborType, 
                       CaseControl=CaseControl, Original = Original, Approximate = Approximate, 
                       Adjust=Adjust, MaxRange=MaxRange, 
                       CheckArguments = FALSE,
                       simulate=SimulatedPP, savefuns=TRUE
                      )
  attr(Envelope, "einfo")$H0 <- switch (SimulationType,
                                        RandomLocation = "Random Location",
                                        RandomLabeling = "Random Labeling",
                                        PopulationIndependence = "Population Independence"
  )
  # Calculate confidence intervals
  Envelope <- FillEnvelope(Envelope, Alpha, Global)
  # No edge effect correction
  attr(Envelope, "einfo")$valname <- NULL
  # Return the envelope
  return (Envelope)
}
