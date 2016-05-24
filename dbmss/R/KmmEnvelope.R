KmmEnvelope <-
function(X, r = NULL, NumberOfSimulations = 100, Alpha = 0.05, ReferenceType = "", Global = FALSE) {

  CheckdbmssArguments()
  
  # The only null hypothesis is random labelling (equivalently, random location)
  SimulatedPP <- expression(rRandomLocation(X, ReferenceType, CheckArguments = FALSE))
  
  # local envelope, keep extreme values for lo and hi (nrank=1)
  Envelope <- envelope(X, fun=Kmmhat, nsim=NumberOfSimulations, nrank=1,
                       r=r, ReferenceType=ReferenceType, 
                       CheckArguments = FALSE,
                       simulate=SimulatedPP, savefuns=TRUE
                       )
  attr(Envelope, "einfo")$H0 <- "Random Location"
  
  # Calculate confidence intervals
  Envelope <- FillEnvelope(Envelope, Alpha, Global)
  # Return the envelope
  return (Envelope)
}
