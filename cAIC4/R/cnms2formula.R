cnms2formula <-
function(cnms) {
  # A function that builds a random effects formula from the “component names”, 
  # a list that can be extracted from an lmerMod object by .@cnms or 
  # getME(., "cnms").
  #
  # Args: 
  #   cnms     = List from an lmerMod object by .@cnms or getME(., "cnms").
  #
  # Returns:
  #   reFormula = random effects part of a lmerMod formula
  #
  len      <- unlist(lapply(cnms, length))
  cnms     <- cnms[which(len != 0)]
  charForm <- character(length(cnms))
  
  for(i in 1:length(cnms)) {
    if (cnms[[i]][1] == "(Intercept)") {
      cnms[[i]][1] <- "1"
    } else {
      cnms[[i]][1] <- "-1"
    }
    charForm[i] <- paste("(", paste(cnms[[i]], collapse = " + "), 
                         " | ",names(cnms)[i], ")", sep = "")
  }
  
  reFormula <- paste(charForm, collapse = " + ")

  return(reFormula)
}
