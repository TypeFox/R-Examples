# Modified: 5 Nov 2015

molMass <- function(form) {

  # Check argument
  checkArgClassValue(form, 'character')

  # Loop through all elements in form
  mmass <- NULL
  for(f in form) {
    # If and only if first letter of form is lowercase, entire string is capitalized
    if(grepl('^[a-z]', f)) f <- toupper(f) 

    # Get coefficients of formula
    fc <- readFormula(f)

    # Check for unidentified element
    if(any(!names(fc) %in% names(atom.weights))) stop('One or more elements in \"form\" is not in the database. You can add it to the \"atom.weights\" vector if you want to modify the function code. Otherwise send a request to saha@kbm.sdu.dk.')

    # Calculate molar mass, using names of fc for indexing
    mmass <- c(mmass, sum(atom.weights[names(fc)]*fc))
  }

  #names(mmass) <- form

  return(mmass)
}
