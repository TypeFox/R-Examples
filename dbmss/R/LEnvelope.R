LEnvelope <-
function(X, r = NULL, NumberOfSimulations = 100, Alpha = 0.05, 
         ReferenceType = "", NeighborType = "", SimulationType = "RandomPosition", Global = FALSE) {
  
  # Calculate the envelope of K
  Envelope <- KEnvelope(X, r, NumberOfSimulations, Alpha, ReferenceType, NeighborType, SimulationType, Global)
  # Transform K to L
  Columns <- names(Envelope)[-1]
  for(i in Columns) {
    Envelope[[i]] <- sqrt(Envelope[[i]]/pi)-Envelope$r
  }
  attr(Envelope, "ylab") <- "L(r)"
  attr(Envelope, "yexp") <- "L(r)"
  attr(Envelope, "fname") <- "L"
  return (Envelope)
}
