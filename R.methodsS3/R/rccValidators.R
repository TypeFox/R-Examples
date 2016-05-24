rccValidateFunctionName <- function(name, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate 'name'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert that the generic function name is a valid function name.
  firstLetter <- substring(gsub("^[.]*", "", name), 1,1);

  allowedFirst <- c("?", "$", "$<-", "[", "[<-", "[[", "[[<-");
  allowedFirst <- c(allowedFirst, "+", "-", "*", "^", "%");
  if (!is.element(firstLetter, allowedFirst)) {
    if (!is.element(tolower(firstLetter), letters))
      throw("Except for a few operators, method/function names must begin with a letter: ", name);

    # Check first letter  
    if (firstLetter == toupper(firstLetter))
      throw("Method/function names should start with a lower case letter: ", name);
  }
}
export(rccValidateFunctionName) <- FALSE;

rccValidateSetMethodS3 <- function(name, ...) {
  rccValidateFunctionName(name=name)
}
export(rccValidateSetMethodS3) <- FALSE;

rccValidateSetGenericS3 <- function(name, ...) {
  rccValidateFunctionName(name=name)
}
export(rccValidateSetGenericS3) <- FALSE;


############################################################################
# HISTORY:
# 2012-06-22
# o Now rccValidateFunctionName() also accepts names starting with
#   symbols "+", "-", "*", "^", and "%".
# 200x-xx-xx
# o Created.
############################################################################
