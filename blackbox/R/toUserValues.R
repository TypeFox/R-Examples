toUserValues <- function(innervalues) {
  innervalues[names(innervalues)=="latt2Ns2"] <- innervalues*blackbox.getOption("Nbfactor")
  innervalues[names(innervalues)=="condS2"] <- innervalues*blackbox.getOption("S2factor")
  return(as.numeric(innervalues)) ## makes sure that the names are lost (since they no longer are all valid: force user to care about this)
}
