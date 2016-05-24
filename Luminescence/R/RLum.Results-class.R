#' @include get_RLum.R set_RLum.R length_RLum.R names_RLum.R
NULL

#' Class \code{"RLum.Results"}
#'
#' Object class contains results data from functions (e.g. \code{\link{analyse_SAR.CWOSL}}.
#'
#' @name RLum.Results-class
#'
#' @docType class
#'
#' @slot data Object of class "list" containing output data
#'
#' @note The class is intended to store results from functions to be used by
#' other functions. The data in the object should always be accessed by the
#' method \code{get_RLum}.
#'
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("RLum.Results", ...)}.
#'
#' @section Class version: 0.4.0
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#'
#' @seealso \code{\linkS4class{RLum}}, \code{\link{plot_RLum}}, \code{\link{merge_RLum}}
#'
#' @keywords classes methods
#'
#' @examples
#'
#' showClass("RLum.Results")
#'
#' ##create an empty object from this class
#' set_RLum(class = "RLum.Results")
#'
#' ##use another function to show how it works
#'
#' ##Basic calculation of the dose rate for a specific date
#'  dose.rate <-  calc_SourceDoseRate(
#'    measurement.date = "2012-01-27",
#'    calib.date = "2014-12-19",
#'    calib.dose.rate = 0.0438,
#'    calib.error = 0.0019)
#'
#' ##show object
#' dose.rate
#'
#' ##get results
#' get_RLum(dose.rate)
#'
#' ##get parameters used for the calcualtion from the same object
#' get_RLum(dose.rate, data.object = "parameters")
#'
#' ##alternatively objects can be accessed using S3 generics, such as
#' dose.rate$parameters
#'
#' @export
setClass(
  Class = "RLum.Results",
  slots = list(data = "list"),
  contains = "RLum",
  prototype = list (data = list())
)


####################################################################################################
###as()
####################################################################################################
##LIST
##COERCE RLum.Results >> list AND list >> RLum.Results
#' as() - RLum-object coercion
#'
#' for \code{[RLum.Results]}
#'
#' \bold{[RLum.Results]}\cr
#'
#' \tabular{ll}{
#'  \bold{from} \tab \bold{to}\cr
#'   \code{list} \tab \code{list}\cr
#' }
#'
#' Given that the \code{\link{list}} consits of \code{\linkS4class{RLum.Results}} objects.
#'
#' @name as
#'
#'
setAs("list", "RLum.Results",
      function(from,to){

        new(to,
            orginator = "coercion",
            data = from)

      })

setAs("RLum.Results", "list",
      function(from){

        from@data

      })

####################################################################################################
###show()
####################################################################################################
#' @describeIn RLum.Results
#' Show structure of \code{RLum.Results} object
#' @export
setMethod("show",
          signature(object = "RLum.Results"),
          function(object){


            ##data elements
            temp.names <- names(object@data)

            if(length(object) > 0){
              temp.type <- sapply(1:length(object@data),
                                  function(x){

                                    paste("\t .. $", temp.names[x],
                                          " : ",
                                          is(object@data[[x]])[1],
                                          sep = "")


                                  })
            }else{
              temp.type <- paste0("\t .. $", temp.names, " : ",is(object@data)[1])

            }

            temp.type <- paste(temp.type, collapse="\n")

            ##print information
            cat("\n [RLum.Results]")
            cat("\n\t originator: ", object@originator,"()", sep="")
            cat("\n\t data:", length(object@data))
            cat("\n", temp.type)


          }
)



####################################################################################################
###set_RLum()
####################################################################################################
#' @describeIn RLum.Results
#' Construction method for an RLum.Results object.
#'
#' @param class [\code{set_RLum}] \code{\link{character}} \bold{(required)}: name of the \code{RLum} class to create
#' @param originator [\code{set_RLum}] \code{\link{character}} (automatic): contains the
#' name of the calling function
#' (the function that produces this object); can be set manually.
#' @param data [\code{set_RLum}] \code{\link{list}} (optional): a list containing the data to
#' be stored in the object
#'
#' @return
#'
#' \bold{\code{set_RLum}}:\cr
#'
#' Returns an object from the class \code{\linkS4class{RLum.Results}}\cr
#'
#' @export
setMethod("set_RLum",
          signature = signature("RLum.Results"),

          function(class, originator, data = list()){

            new(
              Class = "RLum.Results",
              originator = originator,
              data = data)
          })


####################################################################################################
###get_RLum()
####################################################################################################
#' @describeIn RLum.Results
#' Accessor method for RLum.Results object. The argument data.object allows
#' directly accessing objects delivered within the slot data. The default
#' return object depends on the object originator (e.g., \code{fit_LMCurve}).
#' If nothing is specified always the first \code{data.object} will be returned.
#'
#' Note: Detailed specification should be made in combination with the originator slot in the
#' receiving function if results are pipped.
#'
#' @param object [\code{get_RLum}] \code{\linkS4class{RLum.Results}} (required): an object of class
#' \code{\linkS4class{RLum.Results}} to be evaluated
#'
#' @param data.object [\code{get_RLum}] \code{\link{character}} or
#' \code{\link{numeric}}: name or index of the data slot to be returned
#'
#' @param drop [\code{get_RLum}] \code{\link{logical}} (with default): coerce to the next possible layer
#' (which are data objects, \code{drop = FALSE} keeps the original \code{RLum.Results}
#'
#' @return
#'
#' \bold{\code{get_RLum}}:\cr
#'
#' Returns: \cr
#' (1) Data object from the specified slot \cr
#' (2) \code{\link{list}} of data objects from the slots if 'data.object' is vector or \cr
#' (3) an \code{\linkS4class{RLum.Results}} for \code{drop = FALSE}.\cr
#'
#'
#' @export
setMethod("get_RLum",
          signature = signature("RLum.Results"),
          definition = function(object, data.object, drop = TRUE) {

            if (!missing(data.object)) {

              ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              ##CASE1: data.object is of type 'character'
              if (is(data.object, "character")) {

                #check if the provided names are available
                if (all(data.object %in% names(object@data))) {

                  ##account for multiple inputs
                  if (length(data.object) > 1) {
                    temp.return <- sapply(data.object, function(x){
                      object@data[[x]]

                    })

                  } else{
                    temp.return <- list(data.object = object@data[[data.object]])

                  }


                } else{
                  error.message <- paste0(
                    "[get_RLum()] data.object(s) unknown, valid names are: ",
                    paste(names(object@data), collapse = ", ")

                  )
                  stop(error.message)

                }

              }

              ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              ##CASE2: data.object is of type 'numeric'
              else if (is(data.object, "numeric")) {
                ##check if index is valid
                if (max(data.object) > length(object@data)) {
                  stop("[get_RLum] 'data.object' index out of bounds!")

                } else if (length(data.object) > 1) {
                  temp.return <- lapply(data.object, function(x) {
                    object@data[[x]]

                  })


                } else{
                  temp.return <- list(object@data[[data.object]])

                }

                ##restore names as that get los with this method
                names(temp.return) <- names(object@data)[data.object]

              }
              ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              ##CASE3: data.object is of an unsupported type
              else{
                stop("[get_RLum] 'data.object' has to be of type character or numeric!")
              }

            ##the CASE data.object is missing
            } else{

              ##return always the first object if nothing is specified
              temp.return <- object@data[1]

            }

          ##CHECK whether an RLum.Results object needs to be produced ...
          ##This will just be the case if the funtion havn't returned something before
          if (drop) {
            ##we need to access the list here, otherwise we get unexpected behaviour as drop = TRUE
            ##should always return the lowest possible element here
            return(temp.return[[1]])

          } else{

            return(set_RLum(
              "RLum.Results",
              originator = object@originator,
              data = temp.return
            ))


          }
    })



####################################################################################################
###length_RLum()
####################################################################################################
#' @describeIn RLum.Results
#' Returns the length of the object, i.e., number of stored data.objects
#'
#' @return
#'
#' \bold{\code{length_RLum}}\cr
#'
#' Returns the number of data elements in the \code{RLum.Results} object.
#'
#' @export
setMethod("length_RLum",
          "RLum.Results",
          function(object){

            length(object@data)

          })

####################################################################################################
###names_RLum()
####################################################################################################
#' @describeIn RLum.Results
#' Returns the names data.objects
#'
#' @return
#'
#' \bold{\code{names_RLum}}\cr
#'
#' Returns the names of the data elements in the object.
#'
#' @export
setMethod("names_RLum",
          "RLum.Results",
          function(object){
             names(object@data)

          })
