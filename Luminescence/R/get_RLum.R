#' General accessor function for RLum S4 class objects
#'
#' Function calls object-specific get functions for RLum S4 class objects.
#'
#' The function provides a generalised access point for specific
#' \code{\linkS4class{RLum}} objects.\cr Depending on the input object, the
#' corresponding get function will be selected. Allowed arguments can be found
#' in the documentations of the corresponding \code{\linkS4class{RLum}} class.
#'
#' @param object \code{\linkS4class{RLum}} (\bold{required}): S4 object of
#' class \code{RLum} or an object of type \code{\link{list}} containing only objects of type
#' \code{\linkS4class{RLum}}
#'
#' @param \dots further arguments that will be passed to the object specific methods. For
#' furter details on the supported arguments please see the class
#' documentation: \code{\linkS4class{RLum.Data.Curve}},
#' \code{\linkS4class{RLum.Data.Spectrum}}, \code{\linkS4class{RLum.Data.Image}},
#' \code{\linkS4class{RLum.Analysis}} and \code{\linkS4class{RLum.Results}}
#'
#' @return Return is the same as input objects as provided in the list.
#'
#' @section Function version: 0.3.0
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#'
#' @seealso
#' \code{\linkS4class{RLum.Data.Curve}},
#' \code{\linkS4class{RLum.Data.Image}},
#' \code{\linkS4class{RLum.Data.Spectrum}},
#' \code{\linkS4class{RLum.Analysis}},
#' \code{\linkS4class{RLum.Results}}
#'
#' @keywords utilities
#'
#' @aliases get_RLum.Data.Curve get_RLum.Data.Image get_RLum.Data.Spectrum
#' get_RLum.Analysis get_RLum.Results
#'
#' @examples
#'
#'
#' ##Example based using data and from the calc_CentralDose() function
#'
#' ##load example data
#' data(ExampleData.DeValues, envir = environment())
#'
#' ##apply the central dose model 1st time
#' temp1 <- calc_CentralDose(ExampleData.DeValues$CA1)
#'
#' ##get results and store them in a new object
#' temp.get <- get_RLum(object = temp1)
#'
#'
#' @export
setGeneric("get_RLum", function (object, ...) { standardGeneric("get_RLum") })

# Method for get_RLum method for RLum objects in a list for a list of objects  -------------------
#' @describeIn get_RLum
#' Returns a list of \code{\linkS4class{RLum}} objects that had been passed to \code{\link{get_RLum}}
#'
#' @param null.rm \code{\link{logical}} (with default): option to get rid of empty and NULL objects
#'
#' @export
setMethod("get_RLum",
          signature = "list",
          function(object, null.rm = FALSE, ...){

            selection <- lapply(1:length(object), function(x){

              ##get rid of all objects that are not of type RLum, this is better than leaving that
              ##to the user
              if(inherits(object[[x]], what = "RLum")){

                get_RLum(object[[x]],...)

              }else{

               warning(paste0("[get_RLum()] object #",x," in the list was not of type 'RLum' and has been removed!"),
                       call. = FALSE)
                return(NULL)

                return(NULL)
              }

            })

            ##remove empty or NULL objects
            if(null.rm){

                ##first set all empty objects to NULL ... for RLum.Analysis objects
                selection <- lapply(1:length(selection), function(x){
                  if(is(selection[[x]], "RLum.Analysis") && length(selection[[x]]@records) == 0){
                    return(NULL)

                  }else{
                    return(selection[[x]])

                  }

                })

                ##get rid of all NULL objects
                selection <- selection[!sapply(selection, is.null)]


            }

            return(selection)

          })


## ---- DEPRECATED GENERICS
# .Deprecated in package version 0.5.0
# .Defunct in 0.5.1
# Removed in 0.6.0

#' @noRd
#' @export
get_RLum.Analysis <- function(...) {
  .Defunct("get_RLum")
  get_RLum(...)
}

#' @noRd
#' @export
get_RLum.Data.Curve <- function(...) {
  .Defunct("get_RLum")
  get_RLum(...)
}

#' @noRd
#' @export
get_RLum.Data.Image <- function(...) {
  .Defunct("get_RLum")
  get_RLum(...)
}

#' @noRd
#' @export
get_RLum.Data.Spectrum <- function(...) {
  .Defunct("get_RLum")
  get_RLum(...)
}

#' @noRd
#' @export
get_RLum.Results <- function(...) {
  .Defunct("get_RLum")
  get_RLum(...)
}
