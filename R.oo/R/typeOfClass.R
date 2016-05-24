###########################################################################/**
# @RdocDefault typeOfClass
#
# @title "Gets the type of a class (S3 or S4)"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{object}{The object to be checks.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string \code{"S3"}, \code{"S3-Object"} or \code{"S4"},
#  or @NA if neither.
# }
#
# @author
#
# \keyword{character}
#*/###########################################################################
setMethodS3("typeOfClass", "default", function(object, ...) {
  if (is.null(object))
    return(NA_character_);

  if (inherits(object, "classRepresentation"))
    return("S4");

  if (is.character(object)) {
    if (isClass(object)) {
      if (inherits(object, "oldClass"))
        return("S3");
      return("S4");
    }

    # TO DO/FIX ME: This part only works when packages are attached.
    # /HB 2013-10-08
    if (!exists(object, mode="function"))
      return(NA_character_);
    object <- get(object, mode="function");
  }

  if (is.function(object) && inherits(object, "Class"))
    return("S3-Object");

  return(NA_character_);
})

############################################################################
# 2013-10-08
# o Now returning NA_character_ instead of NA.
# 2005-02-15
# o Added arguments '...' in order to match any generic functions.
############################################################################
