#' @include get_RLum.R set_RLum.R names_RLum.R
NULL

#' Class \code{"RLum.Data.Image"}
#'
#' Class for representing luminescence image data (TL/OSL/RF). Such data are for example produced
#' by the function \code{\link{read_SPE2R}}
#'
#' @name RLum.Data.Image-class
#'
#' @docType class
#'
#' @slot recordType Object of class \code{\link{character}}
#' containing the type of the curve (e.g. "OSL image", "TL image")
#'
#' @slot curveType Object of class \code{\link{character}} containing curve type, allowed values
#' are measured or predefined
#'
#' @slot data Object of class \code{\link[raster]{brick}} containing images (raster data).
#'
#' @slot info Object of class \code{\link{list}} containing further meta information objects
#'
#' @note The class should only contain data for a set of images. For additional
#' elements the slot \code{info} can be used.
#'
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{set_RLum("RLum.Data.Image", ...)}.
#'
#' @section Class version: 0.3.0
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne (France)
#'
#' @seealso \code{\linkS4class{RLum}}, \code{\linkS4class{RLum.Data}},
#' \code{\link{plot_RLum}}, \code{\link{read_SPE2R}}
#'
#' @keywords classes
#'
#' @examples
#'
#' showClass("RLum.Data.Image")
#'
#' ##create empty RLum.Data.Image object
#' set_RLum(class = "RLum.Data.Image")
#'
#' @importClassesFrom raster RasterBrick
#' @export
setClass(
  "RLum.Data.Image",
  slots = list(
    recordType = "character",
    curveType = "character",
    data = "RasterBrick",
    info = "list"
  ),
  contains = "RLum.Data",
  prototype = list (
    recordType = character(),
    curveType = character(),
    data = raster::brick(raster::raster(matrix())),
    info = list()
  )
)


####################################################################################################
###as()
####################################################################################################

##DATA.FRAME
##COERCE RLum.Data.Image >> data.frame AND data.frame >> RLum.Data.Image
#' as()
#'
#' for \code{[RLum.Data.Image]}
#'
#' \bold{[RLum.Data.Image]}\cr
#'
#' \tabular{ll}{
#'  \bold{from} \tab \bold{to}\cr
#'   \code{data.frame} \tab \code{data.frame}\cr
#'   \code{matrix} \tab \code{matrix}
#'
#' }
#'
#' @name as
#'
#'
setAs("data.frame", "RLum.Data.Image",
      function(from,to){

        new(to,
            recordType = "unkown curve type",
            curveType = "NA",
            data = as.matrix(from),
            info = list())
      })

setAs("RLum.Data.Image", "data.frame",
      function(from){

        data.frame(x = from@data@values[seq(1,length(from@data@values), by = 2)],
                   y = from@data@values[seq(2,length(from@data@values), by = 2)])

      })


##MATRIX
##COERCE RLum.Data.Image >> matrix AND matrix >> RLum.Data.Image
setAs("matrix", "RLum.Data.Image",
      function(from,to){

        new(to,
            recordType = "unkown curve type",
            curveType = "NA",
            data = raster::brick(raster::raster(as.matrix(from))),
            info = list())
      })

setAs("RLum.Data.Image", "matrix",
      function(from){

        ##only the first object is convertec
        as.matrix(from[[1]])

      })


####################################################################################################
###show()
####################################################################################################
#' @describeIn RLum.Data.Image
#' Show structure of \code{RLum.Data.Image} object
#' @export
setMethod("show",
          signature(object = "RLum.Data.Image"),
          function(object){

            x.rows <- object@data@ncols
            y.cols <- object@data@nrows
            z.range <- paste(min(object@data@data@min),":",max(object@data@data@max))

            ##print information

            cat("\n [RLum.Data.Image]")
            cat("\n\t recordType:", object@recordType)
            cat("\n\t curveType:",  object@curveType)
            cat("\n\t .. recorded frames:", length(object@data@data@names))
            cat("\n\t .. .. pixel per frame:", x.rows*y.cols)
            cat("\n\t .. .. x dimension [px]:", x.rows)
            cat("\n\t .. .. y dimension [px]:", y.cols)
            cat("\n\t .. .. full pixel value range:", z.range)
            cat("\n\t additional info elements:", length(object@info))
            #cat("\n\t\t >> names:", names(object@info))
          }
)


####################################################################################################
###set_RLum()
####################################################################################################
#' @describeIn RLum.Data.Image
#' Construction method for RLum.Data.Image object. The slot info is optional
#' and predefined as empty list by default..
#'
#' @param class \code{[set_RLum]}\code{\link{character}}: name of the \code{RLum} class to create
#' @param originator \code{[set_RLum]} \code{\link{character}} (automatic):
#' contains the name of the calling function (the function that produces this object); can be set manually.
#' @param recordType \code{[set_RLum]} \code{\link{character}}: record type (e.g. "OSL")
#' @param curveType \code{[set_RLum]} \code{\link{character}}: curve type (e.g. "predefined" or "measured")
#' @param data \code{[set_RLum]} \code{\link{matrix}}: raw curve data. If data is of type \code{RLum.Data.Image}
#' this can be used to re-construct the object.
#' @param info \code{[set_RLum]} \code{\link{list}}: info elements
#'
#' @return
#'
#' \bold{\code{set_RLum}}\cr
#'
#' Returns an object from class \code{RLum.Data.Image}
#'
#' @export
setMethod("set_RLum",
          signature = signature("RLum.Data.Image"),

          definition = function(class,
                                originator,
                                recordType = "Image",
                                curveType = NA_character_,
                                data = raster::brick(raster::raster(matrix())),
                                info = list()){

            if (is(data, "RLum.Data.Image")) {

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
              Class = "RLum.Data.Image",
              originator = originator,
              recordType = recordType,
              curveType = curveType,
              data = data@data,
              info = info
            )

            }else{

              new(
                Class = "RLum.Data.Image",
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
#' @describeIn RLum.Data.Image
#' Accessor method for RLum.Data.Image object. The argument info.object is
#'  optional to directly access the info elements. If no info element name is
#'  provided, the raw image data (RasterBrick) will be returned.
#'
#' @param object \code{[show_RLum]}\code{[get_RLum]}\code{[names_RLum]} an object
#' of class \code{\linkS4class{RLum.Data.Image}}
#' @param info.object \code{[get_RLum]} \code{\link{character}} name of the info object to returned
#'
#' @return
#'
#' \bold{\code{get_RLum}}\cr
#'
#' (1) Returns the data object (\code{\link[raster]{brick}})\cr
#' (2) only the info object if \code{info.object} was set.\cr
#'
#' @export
setMethod("get_RLum",
          signature("RLum.Data.Image"),
          definition = function(object, info.object) {

            ##Check if function is of type RLum.Data.Image
            if(is(object, "RLum.Data.Image") == FALSE){

              stop("[get_RLum] Function valid for 'RLum.Data.Image' objects only!")

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
#' @describeIn RLum.Data.Image
#' Returns the names info elements coming along with this curve object
#'
#' @return
#'
#' \bold{\code{names_RLum}}\cr
#'
#' Returns the names of the info elements
#'
#' @export
setMethod("names_RLum",
          "RLum.Data.Image",
          function(object) {
            names(object@info)

          })
