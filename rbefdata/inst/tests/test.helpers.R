# This can hold helpers maybe also some sort of factories to create a certain
# condition for testing setting the right options. This helper needs to be
# sourced into test files that make use of it with:

# source("test.helpers.R")

given_the_user_is <- function(condition = NULL) {
  if(condition == "valid") {
    valid_credentials = "i9ZtE0F5qoma1EooxVM"
    bef.options(user_credentials = valid_credentials)
    return(valid_credentials)
  }
  if(condition == "invalid") {
    bef.options(user_credentials="invalid")
    return("invalid")
  } else {
    warning("Valid parameters are: valid and invalid")
  }
}

given_the_portal_is <- function(environment = NULL) {
  if(environment == "development") {
    bef.options(url = "http://befdatadevelopment.biow.uni-leipzig.de")
    return("development")
  }
  if(environment == "production") {
    bef.options(url = "http://befdataproduction.biow.uni-leipzig.de")
    return("production")
  } else {
    warning("Valid parameters are: development and production")
  }
}

given_the_dataset_is <- function(available = NULL) {
  if(available == TRUE) {
    id = 1
    return(id)
  }
  if(available == FALSE) {
    id = 11230982734
    return(id)
  } else {
    warning("Valid parameters are: TRUE and FALSE")
  }
}


