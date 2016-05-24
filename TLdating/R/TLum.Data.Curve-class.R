#' Class \code{"TLum.Data.Curve"}
#'
#' Class for luminescence curve data.
#'
#'
#' @name TLum.Data.Curve-class
#' @rdname TLum.Data.Curve-class
#'
#' @slot recordType
#'  \link{character}: record type.
#' @slot curveType
#'  \link{character}: curve type
#' @slot metadata
#'  \link{list}: metadata elements.
#' @slot temperatures
#'  \link{numeric}: Object containing numeric vector with temperature data.
#' @slot data
#'  \link{numeric}: Object containing numeric vector with count data
#' @slot error
#'  \link{numeric}: Object containing numeric vector with count data absolute uncertainty.
#' @slot analysis
#'  \link{list}: data produced by analysis function.
#' @slot .RESERVED
#'  \link{list}: Object containing list of undocumented raw values for internal use only.
#'
#' @note The code and the structure of this class is based on the \linkS4class{RLum.Data.Curve} class from the \link{Luminescence} package.
#'
#' @keywords classes
#'
#' @author David Strebler
#'
#' @exportClass TLum.Data.Curve

setClass(Class = "TLum.Data.Curve",
         contains = "TLum.Data",
         slots = c(recordType = "character",
                   curveType = "character",
                   metadata = "list",
                   temperatures = "numeric",
                   data = "numeric",
                   error= "numeric",
                   analysis = "list",
                   .RESERVED = "list"),
         prototype = list(recordType = character(),
                          curveType = character(),
                          metadata = list(),
                          temperatures = vector(mode="numeric"),
                          data = vector(mode="numeric"),
                          error = vector(mode="numeric"),
                          analysis = list(),
                          .RESERVED = list()
                          )
         )

# Show

#' @rdname TLum.Data.Curve-class
#' @aliases show,TLum.Data.Curve-method

setMethod(f = "show",
          signature = "TLum.Data.Curve",
          definition = function(object){

            cat("\n [TLum.Data.Curve]")
            cat("\n\t recordType:", object@recordType)
            cat("\n\t curveType:",  object@curveType)
            cat("\n\t measured values:", length(object@data))
            cat("\n\t .. range of temperatures:", range(object@temperatures))
            cat("\n\t .. range of values:", range(object@data))
            cat("\n\t .. range of uncertainties:", range(object@error))
            cat("\n\t additional information:", length(object@metadata))
            cat("\n\t additional analysis data:", length(object@analysis))

          })

# set

#' @name TLum.Data.Curve-class
#' @rdname TLum.Data.Curve-class
#'
#' @param recordType
#'  \link{character}: record type.
#' @param curveType
#'  \link{character}: curve type
#' @param metadata
#'  \link{list}: metadata elements.
#' @param temperatures
#'  \link{numeric}: Object containing numeric vector with temperature data.
#' @param data
#'  \link{numeric}: Object containing numeric vector with count data
#' @param error
#'  \link{numeric}: Object containing numeric vector with count data absolute uncertainty.
#' @param analysis
#'  \link{list}: data produced by analysis function.
#' @param .RESERVED
#'  \link{list}: Object containing list of undocumented raw values for internal use only.
#'
#' @exportMethod set_TLum.Data.Curve
setGeneric(name = "set_TLum.Data.Curve",
           def = function(recordType, curveType, metadata, temperatures, data, error, analysis, .RESERVED) {standardGeneric("set_TLum.Data.Curve")}
           )

#' @rdname TLum.Data.Curve-class
#' @aliases set_TLum.Data.Curve set_TLum.Data.Curve,TLum.Data.Curve-method

setMethod(f="set_TLum.Data.Curve",
          signature = c(recordType = "ANY",
                        curveType = "ANY",
                        metadata = "ANY",
                        temperatures= "ANY",
                        data = "ANY",
                        error = "ANY",
                        analysis = "ANY",
                        .RESERVED = "ANY"
                        ),
          definition = function(recordType, curveType, metadata, temperatures, data, error, analysis,  .RESERVED){

            if(missing(recordType)){
              stop("[set_TLum.Data.Curve] Error: Input 'recordType' is missing.")
            }else if(!is.character(recordType)){
              stop("[set_TLum.Data.Curve] Error: Input 'recordType' is not of type 'character'.")
            }

            if(missing(curveType)){
              curveType <- "NA"
            }else if(!is.character(curveType)){
              stop("[set_TLum.Data.Curve] Error: Input 'curveType' is not of type 'character'.")
            }

            if(missing(temperatures)){
              stop("[set_TLum.Data.Curve] Error: Input 'temperatures' is missing.")
            }else if(!is.numeric(temperatures)){
              stop("[set_TLum.Data.Curve] Error: Input 'temperatures' is not of type 'numeric'.")
            }

            if(missing(data)){
              stop("[set_TLum.Data.Curve] Error: Input 'data' is missing.")
            }else if(!is.numeric(data)){
              stop("[set_TLum.Data.Curve] Error: Input 'data' is not of type 'numeric'.")
            }

            if(missing(error)){
              stop("[set_TLum.Data.Curve] Error: Input 'error' is missing.")
            }else if(!is.numeric(error)){
              stop("[set_TLum.Data.Curve] Error: Input 'error' is not of type 'numeric'.")
            }

            if(missing(metadata)){
              metadata <- list()
            }else if(!is.list(metadata)){
              stop("[set_TLum.Data.Curve] Error: Input 'metadata' is not of type 'list'.")
            }
            if(length(metadata) < 1 ){
              warning("[set_TLum.Data.Curve] Warning: Input 'metadata' is empty.")
            }

            if(missing(analysis)){
              analysis <- list()
            }else if(!is.list(analysis)){
              stop("[set_TLum.Data.Curve] Error: Input 'analysis' is not of type 'list'.")
            }

            if(missing(.RESERVED)){
              .RESERVED <- list()
            }else if(!is.list(.RESERVED)){
              stop("[set_TLum.Data.Curve] Error: Input '.RESERVED' is not of type 'list'.")
            }

            new(Class = "TLum.Data.Curve",
                recordType=recordType,
                curveType=curveType,
                metadata=metadata,
                temperatures=temperatures,
                data=data,
                error=error,
                analysis=analysis,
                .RESERVED = .RESERVED)
          })

#get

#' @name TLum.Data.Curve-class
#' @rdname TLum.Data.Curve-class
#'
#' @param object
#'  \linkS4class{TLum.Data.Curve}:  an object of class 'TLum.Data.Curve'.
#' @param ref
#'  \link{character}: name of the wanted element.
#'
#' @exportMethod get_TLum.Data.Curve
setGeneric(name = "get_TLum.Data.Curve",
           def = function(object, ref) {standardGeneric("get_TLum.Data.Curve")}
           )

#' @rdname TLum.Data.Curve-class
#' @aliases get_TLum.Data.Curve get_TLum.Data.Curve,TLum.Data.Curve-method

setMethod(f = "get_TLum.Data.Curve",
          signature = c(object="ANY",
                        ref = "ANY"),
          definition = function(object, ref){

            if(!is(object, "TLum.Data.Curve")){
              stop("[get_TLum.Data.Curve] Error: Function valids for 'TLum.Data.Curve' objects only!")
            }

            if(missing(ref)){
              res <- object@data

            }else{
              if(!is.character(ref)){
                stop("[get_TLum.Data.Curve] Error: Input 'ref' is not of type 'character'.")
              }

              if(ref %in% names(object@metadata)){
                res <- unlist(object@metadata[ref])

              }else if(ref %in% names(object@analysis)){
                res <- unlist(object@analysis[ref])

              }else if(ref == "data" || ref == "values"){
                res <- object@data

              }else if(ref == "error" || ref == "uncertainty"){
                res <- object@error

              }else if(ref == "temperatures"){
                res <- object@temperatures

              }else if(ref == "metadata"){
                res <- object@metadata

              }else if(ref == "analysis"){
                res <- object@analysis

              }else if(ref== ".RESERVED"){
                res <- object@.RESERVED

              }else{
                stop("[get_TLum.Data.Curve] Error: Input 'ref' is unknown.")
              }

              return(res)
            }
          })
