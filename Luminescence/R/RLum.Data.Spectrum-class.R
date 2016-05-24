#' @include get_RLum.R set_RLum.R names_RLum.R
NULL

#' Class \code{"RLum.Data.Spectrum"}
#'
#' Class for representing luminescence spectra data (TL/OSL/RF).
#'
#' @name RLum.Data.Spectrum-class
#'
#' @docType class
#'
#' @slot recordType Object of class \code{\link{character}} containing the type of the curve (e.g. "TL" or "OSL")
#'
#' @slot curveType Object of class \code{\link{character}} containing curve type, allowed values
#' are measured or predefined
#'
#' @slot data Object of class \code{\link{matrix}} containing spectrum (count) values.
#' Row labels indicate wavelength/pixel values, column labels are temperature or time values.
#'
#' @slot info Object of class \code{\link{list}} containing further meta information objects
#'
#' @note The class should only contain data for a single spectra data set. For
#' additional elements the slot \code{info} can be used. Objects from this class are automatically
#' created by, e.g., \code{\link{read_XSYG2R}}
#'
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{set_RLum("RLum.Data.Spectrum", ...)}.
#'
#' @section Class version: 0.3.0
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne (France)
#'
#' @seealso \code{\linkS4class{RLum}}, \code{\linkS4class{RLum.Data}},
#' \code{\link{plot_RLum}}
#'
#' @keywords classes
#'
#' @examples
#'
#' showClass("RLum.Data.Spectrum")
#'
#' ##show example data
#' data(ExampleData.XSYG, envir = environment())
#' TL.Spectrum
#'
#' ##show data matrix
#' get_RLum(TL.Spectrum)
#'
#' ##plot spectrum
#' \dontrun{
#' plot_RLum(TL.Spectrum)
#' }
#' @export
setClass(
  "RLum.Data.Spectrum",
  slots = list(
    recordType = "character",
    curveType = "character",
    data = "matrix",
    info = "list"
  ),
  contains = "RLum.Data",
  prototype = list (
    recordType = NA_character_,
    curveType = NA_character_,
    data = matrix(),
    info = list()
  )
)


####################################################################################################
###as()
####################################################################################################
##data.frame
##COERCE RLum.Data.Spectrum >> data.frame AND data.frame >> RLum.Data.Spectrum
#' as()
#'
#' for \code{[RLum.Data.Spectrum]}
#'
#'
#' \bold{[RLum.Data.Spectrum]}\cr
#'
#' \tabular{ll}{
#'  \bold{from} \tab \bold{to}\cr
#'   \code{data.frame} \tab \code{data.frame}\cr
#'   \code{matrix} \tab \code{matrix}
#'
#' }
#'
#'
#' @name as
#'
#'
setAs("data.frame", "RLum.Data.Spectrum",
      function(from,to){

        new(to,
            recordType = NA_character_,
            curveType = NA_character_,
            data = as.matrix(from),
            info = list())
      })

setAs("RLum.Data.Spectrum", "data.frame",
      function(from){
        as.data.frame(from@data)

      })


##MATRIX
##COERCE RLum.Data.Spectrum >> matrix AND matrix >> RLum.Data.Spectrum
setAs("matrix", "RLum.Data.Spectrum",
      function(from,to){
        new(to,
            recordType = NA_character_,
            curveType = NA_character_,
            data = from,
            info = list())
      })

setAs("RLum.Data.Spectrum", "matrix",
      function(from){
        from@data

      })


####################################################################################################
###show()
####################################################################################################
#' @describeIn RLum.Data.Spectrum
#' Show structure of \code{RLum.Data.Spectrum} object
#' @export
setMethod("show",
          signature(object = "RLum.Data.Spectrum"),
          function(object){

            x.range <- suppressWarnings(range(as.numeric(rownames(object@data))))
            y.range <- suppressWarnings(range(as.numeric(colnames(object@data))))
            z.range <- range(object@data)

            ##print information

            cat("\n [RLum.Data.Spectrum]")
            cat("\n\t recordType:", object@recordType)
            cat("\n\t curveType:",  object@curveType)
            cat("\n\t .. recorded frames:", length(object@data[1,]))
            cat("\n\t .. .. measured values per frame:", length(object@data[,1]))
            cat("\n\t .. .. range wavelength/pixel:", x.range)
            cat("\n\t .. .. range time/temp.:", y.range)
            cat("\n\t .. .. range count values:", z.range)
            cat("\n\t additional info elements:", length(object@info))
            #cat("\n\t\t >> names:", names(object@info))
          }
)


####################################################################################################
###set_RLum()
####################################################################################################
#' @describeIn RLum.Data.Spectrum
#' Construction method for RLum.Data.Spectrum object. The slot info is optional
#'  and predefined as empty list by default
#'
#' @param class [\code{set_RLum}] \code{\link{character}} (automatic): name of the \code{RLum} class to create.
#' @param originator \code{\link{character}} (automatic): contains the name of the calling function
#' (the function that produces this object); can be set manually.
#' @param recordType [\code{set_RLum}] \code{\link{character}}: record type (e.g. "OSL")
#' @param curveType [\code{set_RLum}] \code{\link{character}}: curve type (e.g. "predefined" or "measured")
#' @param data [\code{set_RLum}] \code{\link{matrix}}: raw curve data. If data is of
#' type \code{RLum.Data.Spectrum}, this can be used to re-construct the object.
#' @param info [\code{set_RLum}] \code{\link{list}}: info elements
#'
#' @return
#'
#' \bold{\code{[set_RLum]}}\cr
#'
#' An object from the class \code{RLum.Data.Spectrum}
#'
#' @export
setMethod("set_RLum",
          signature = signature("RLum.Data.Spectrum"),
          definition = function(
            class,
            originator,
            recordType = "Spectrum",
            curveType = NA_character_,
            data = matrix(),
            info = list()){

            if (is(data, "RLum.Data.Spectrum")) {

              ##check for missing curveType
              if (missing(curveType)) {
                curveType <- data@curveType

              }

              ##check for missing recordType
              if(missing(recordType)){
                recordType <- data@recordType

              }

              ##check for missing data ... not possible as data is the object itself

              ##check for missing info
              if(missing(info)){
                info <- data@info

              }


              new(
                Class = "RLum.Data.Spectrum",
                originator = originator,
                recordType = recordType,
                curveType = curveType,
                data = data@data,
                info = info
              )

            }else{
            ##construct object
            new(
              Class = "RLum.Data.Spectrum",
              originator = originator,
              recordType = recordType,
              curveType = curveType,
              data = data,
              info = info
            )

            }
          })

####################################################################################################
###get_RLum()
####################################################################################################
#' @describeIn RLum.Data.Spectrum
#' Accessor method for RLum.Data.Spectrum object. The argument info.object
#' is optional to directly access the info elements. If no info element name
#' is provided, the raw curve data (matrix) will be returned
#'
#' @param object [\code{show_RLum}][\code{get_RLum}][\code{names_RLum}] an object of
#'  class \code{\linkS4class{RLum.Data.Spectrum}}
#' @param info.object [\code{get_RLum}] \code{\link{character}} (optional): the name of the info
#' object to be called
#'
#' @return
#'
#' \bold{\code{get_RLum}}\cr
#'
#' (1) A \code{\link{matrix}} with the spectrum values or \cr
#' (2) only the info object if \code{info.object} was set.\cr
#'
#' @export
setMethod("get_RLum",
          signature("RLum.Data.Spectrum"),
          definition = function(object, info.object) {

            ##Check if function is of type RLum.Data.Spectrum
            if(is(object, "RLum.Data.Spectrum") == FALSE){

              stop("[get_RLum] Function valid for 'RLum.Data.Spectrum' objects only!")

            }

            ##if missing info.object just show the curve values

            if(missing(info.object) == FALSE){

              if(is(info.object, "character") == FALSE){
                stop("[get_RLum] 'info.object' has to be a character!")
              }

              if(info.object %in% names(object@info) == TRUE){

                unlist(object@info[info.object])

              }else{

                ##grep names
                temp.element.names <- paste(names(object@info), collapse = ", ")

                stop.text <- paste("[get_RLum] Invalid element name. Valid names are:", temp.element.names)

                stop(stop.text)

              }


            }else{

              object@data

            }
          })


####################################################################################################
###names_RLum()
####################################################################################################
#' @describeIn RLum.Data.Spectrum
#' Returns the names info elements coming along with this curve object
#'
#' @return
#'
#' \bold{\code{names_RLum}}\cr
#'
#' The names of the info objects
#'
#' @export
setMethod("names_RLum",
          "RLum.Data.Spectrum",
          function(object){
            names(object@info)

          })
