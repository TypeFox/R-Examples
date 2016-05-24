#' @include normEndp.R
NULL

###############################################
# --------------------------------------------#
# Class endpoint                              #
# --------------------------------------------#
###############################################

# All endpoints should be added to this class union

# Common representation of the endpoints.
#
# @name endpoint
setClassUnion("endpoint", "normEndp")


# --------------------------------------------
# Accessor functions for endpoints
# --------------------------------------------

#' Access the expectation value slot of a normEndp S4 object
#' 
#' @param obj object of class normEndp
mu <- function(obj) {
  if (.hasSlot(obj, "mu")) obj@mu else stop("object has no slot named mu.")  
}

#' Function returning the standard deviation slot of a normEndp S4 object
#' 
#' @param obj object of class normEndp
sigma <- function(obj) {
  if (.hasSlot(obj, "sigma")) obj@sigma else stop("object has no slot named sigma.") 
}

#' Method defining the $ operator for the endpoint class
#' 
#' @inheritParams overview
setMethod("$", "endpoint",
          function(x, name) slot(x, name))


# --------------------------------------------
# Show function for endopoints
# --------------------------------------------

setMethod("show", signature = "endpoint", definition = function(object){
  validObject(object)
  # headline
  cat("\n Object of class \"", is(object)[1], "\"\n\n", sep="")
  # iterate through all slots of the object
  names <- slotNames(object)
  for(name in names){
    cat("\t", name, "=", slot(object, name), "\n")
  }
  cat("\n")
})
