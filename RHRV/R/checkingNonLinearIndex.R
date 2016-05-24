############################## checkingNonLinearIndex ##########################
checkingNonLinearIndex <- function(indexNonLinearAnalysis, numberOfIndex){
  if (indexNonLinearAnalysis == -1 ) {
    stop("  --- Nonlinear analysis index not present in function arguments ---\n    --- Quitting now!! ---\n")
  }
  
  if (numberOfIndex < indexNonLinearAnalysis)  {
    stop("   --- NonlinearAnalysis analysis no. ",indexNonLinearAnalysis," not present!! ---\n    --- Quitting now!! ---\n")
  }
  
  if (indexNonLinearAnalysis < 1 )  {
    stop("   --- There are no NonlinearAnalysis structures!!---\n --- Create some using CreateNonLinearAnalysis() ---\n")
  }
}