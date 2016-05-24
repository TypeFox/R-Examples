#' @include get_RLum.R set_RLum.R length_RLum.R structure_RLum.R names_RLum.R
NULL

#' Class \code{"RLum.Analysis"}
#'
#' Object class to represent analysis data for protocol analysis, i.e. all curves, spectra etc.
#' from one measurements. Objects from this class are produced, by e.g. \code{\link{read_XSYG2R}},
#' \code{\link{read_Daybreak2R}}
#'
#'
#' @name RLum.Analysis-class
#'
#' @docType class
#'
#' @slot protocol Object of class \code{\link{character}} describing the applied measurement protocol
#'
#' @slot records Object of class \code{\link{list}} containing objects of class \code{\linkS4class{RLum.Data}}
#'
#' @note The method \code{\link{structure_RLum}} is currently just avaiblable for objects
#' containing \code{\linkS4class{RLum.Data.Curve}}.
#'
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{set_RLum("RLum.Analysis", ...)}.
#'
#' @section Class version: 0.3.0
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#'
#' @seealso \code{\link{Risoe.BINfileData2RLum.Analysis}},
#' \code{\linkS4class{Risoe.BINfileData}}, \code{\linkS4class{RLum}}
#'
#' @keywords classes methods
#'
#' @examples
#'
#' showClass("RLum.Analysis")
#'
#' ##set empty object
#' set_RLum(class = "RLum.Analysis")
#'
#' ###use example data
#' ##load data
#' data(ExampleData.RLum.Analysis, envir = environment())
#'
#' ##show curves in object
#' get_RLum(IRSAR.RF.Data)
#'
#' ##show only the first object, but by keeping the object
#' get_RLum(IRSAR.RF.Data, record.id = 1, drop = FALSE)
#'
#' @export
setClass("RLum.Analysis",
         slots = list(
           protocol = "character",
           records = "list"
         ),
         contains = "RLum",
         prototype = list (
           protocol = NA_character_,
           records = list()
         )
)


####################################################################################################
###as()
####################################################################################################
##LIST
##COERCE RLum.Analyse >> list AND list >> RLum.Analysis
#' as() - RLum-object coercion
#'
#' for \code{[RLum.Analysis]}
#'
#' \bold{[RLum.Analysis]}\cr
#'
#' \tabular{ll}{
#'  \bold{from} \tab \bold{to}\cr
#'   \code{list} \tab \code{list}\cr
#' }
#'
#' Given that the \code{\link{list}} consits of \code{\linkS4class{RLum.Analysis}} objects.
#'
#' @name as
#'
#'
setAs("list", "RLum.Analysis",
      function(from,to){

        new(to,
            protocol = NA_character_,
            records = from)
      })

setAs("RLum.Analysis", "list",
      function(from){

        lapply(1:length(from@records), function(x){
          from@records[[x]]

        })

      })


####################################################################################################
###show()
####################################################################################################
#' @describeIn RLum.Analysis
#' Show structure of \code{RLum.Analysis} object
#' @export
setMethod("show",
          signature(object = "RLum.Analysis"),
          function(object){

            ##print
            cat("\n [RLum.Analysis]")

            ##show slot originator, for compatibily reasons with old example data, here
            ##a check
            if(.hasSlot(object, "originator")){cat("\n\t originator:", paste0(object@originator,"()"))}

            cat("\n\t protocol:", object@protocol)
            cat("\n\t number of records:", length(object@records))

            #skip this part if nothing is included in the object
            if(length(object@records) > 0){

              ##get object class types
              temp <- sapply(1:length(object@records), function(x){

                is(object@records[[x]])[1]

              })

              ##print object class types
              sapply(1:length(table(temp)), function(x){

                ##show RLum class type
                cat("\n\t .. :",names(table(temp)[x]),":",table(temp)[x])


                ##show structure
                ##set width option ... just an implementation for the tutorial output
                ifelse(getOption("width")<=50, temp.width <- 4, temp.width  <- 10)

                cat("\n\t .. .. : ",
                    unlist(sapply(1:length(object@records),  function(i) {

                      if(names(table(temp)[x]) == is(object@records[[i]])[1]){
                        paste(object@records[[i]]@recordType,
                              if(i%%temp.width==0 & i!=length(object@records)){"\n\t .. .. : "})
                      }
                    })))

              })

            }else{

              cat("\n\t >> This is an empty object and cannot be used for further analysis! <<")

            }


          }
)##end show method


####################################################################################################
###set_RLum()
####################################################################################################
#' @describeIn RLum.Analysis
#' Construction method for \code{\linkS4class{RLum.Analysis}} objects.
#'
#' @param class [\code{set_RLum}] \code{\link{character}} (\bold{required}): name of the \code{RLum} class to be created
#' @param originator [\code{set_RLum}] \code{\link{character}} (automatic): contains the name
#' of the calling function (the function that produces this object); can be set manually.
#' @param protocol [\code{set_RLum}] \code{\link{character}} (optional): sets protocol type for
#' analysis object. Value may be used by subsequent analysis functions.
#' @param records [\code{set_RLum}] \code{\link{list}} (\bold{required}): list of \code{\linkS4class{RLum.Analysis}} objects
#'
#' \bold{\code{set_RLum}}:\cr
#'
#' Returns an \code{\linkS4class{RLum.Analysis}} object.
#'
#' @export
setMethod("set_RLum",
          signature = "RLum.Analysis",

          definition = function(class,
                                originator,
                                protocol = NA_character_,
                                records = list()) {


            new(
              Class = "RLum.Analysis",
              originator = originator,
              protocol = protocol,
              records = records
            )


          })

####################################################################################################
###get_RLum()
####################################################################################################
#' @describeIn RLum.Analysis
#' Accessor method for RLum.Analysis object.
#'
#' The slots record.id, recordType, curveType and RLum.type are optional to allow for records
#' limited by their id (list index number), their record type (e.g. recordType = "OSL")
#' or object type.
#'
#' Example: curve type (e.g. curveType = "predefined" or curveType ="measured")
#'
#' The selection of a specific RLum.type object superimposes the default selection.
#' Currently supported objects are: RLum.Data.Curve and RLum.Data.Spectrum
#'
#' @param object \code{[show_RLum]}\code{[get_RLum]}\code{[names_RLum]}\code{[length_RLum]}
#' \code{[structure_RLum]}] an object of class \code{\linkS4class{RLum.Analysis}}
#' (\bold{required})
#'
#' @param record.id [\code{get_RLum}] \code{\link{numeric}} or \code{\link{logical}} (optional): IDs of specific records.
#' If of type \code{logical} the entire id range is assuemd and \code{TRUE} and \code{FALSE} indicates the selection.
#'
#' @param recordType [\code{get_RLum}] \code{\link{character}} (optional): record type (e.g., "OSL").
#' Can be also a vector, for multiple matching, e.g., \code{recordType = c("OSL", "IRSL")}
#'
#' @param curveType [\code{get_RLum}] \code{\link{character}} (optional): curve
#' type (e.g. "predefined" or "measured")
#'
#' @param RLum.type [\code{get_RLum}] \code{\link{character}} (optional): RLum object type.
#' Defaults to "RLum.Data.Curve" and "RLum.Data.Spectrum".
#'
#' @param get.index [\code{get_RLum}] \code{\link{logical}} (optional): return a numeric
#' vector with the index of each element in the RLum.Analysis object.
#'
#' @param recursive [\code{get_RLum}] \code{\link{logical}} (with default): if \code{TRUE} (the default)
#' and the result of the 'get_RLum' request is a single object this object will be unlisted, means
#' only the object itself and no list containing exactly one object is returned. Mostly this makes things
#' easier, however, if this method is used within a loop this might undesired.
#'
#' @param drop [\code{get_RLum}] \code{\link{logical}} (with default): coerce to the next possible layer
#' (which are \code{RLum.Data}-objects), \code{drop = FALSE} keeps the original \code{RLum.Analysis}
#'
#' @return
#'
#' \bold{\code{get_RLum}}:\cr
#'
#' Returns: \cr
#' (1) \code{\link{list}} of \code{\linkS4class{RLum.Data}} objects or \cr
#' (2) Single \code{\linkS4class{RLum.Data}} object, if only one object is contained and
#' \code{recursive = FALSE} or\cr
#' (3) \code{\linkS4class{RLum.Analysis}} ojects for \code{drop = FALSE} \cr
#'
#' @export
setMethod("get_RLum",
          signature = ("RLum.Analysis"),

          function(object, record.id = NULL, recordType = NULL, curveType = NULL, RLum.type = NULL,
                   protocol = "UNKNOWN", get.index = NULL, drop = TRUE, recursive = TRUE){

            ##record.id
            if (is.null(record.id)) {
              record.id <- c(1:length(object@records))

            }else if (!is(record.id, "numeric") & !is(record.id, "logical")) {
              stop("[get_RLum()] 'record.id' has to be of type 'numeric' or 'logical'!")

            }
              ##logical needs a slightly different treatment
              ##Why do we need this? Because a lot of standard R functions work with logical
              ##values instead of numerical indicies
              if(is(record.id, "logical")){
                record.id <- c(1:length(object@records))[record.id]

              }

            ##check if record.id exists
            if (FALSE %in% (abs(record.id) %in% (1:length(object@records)))) {
              stop("[get_RLum()] At least one 'record.id' is invalid!")

            }

            ##recordType
            if (is.null(recordType)) {
              recordType <- unique(unlist(lapply(1:length(object@records),
                                                 function(x) {
                                                   object@records[[x]]@recordType
                                                 })))

            }else{
              if (!is(recordType, "character")) {
                stop("[get_RLum()] 'recordType' has to be of type 'character'!")

              }

            }

            ##curveType
            if(is.null(curveType)) {
              curveType <- unique(unlist(lapply(1:length(object@records),
                                                function(x) {
                                                  object@records[[x]]@curveType
                                                })))

            }else if (!is(curveType, "character")) {
              stop("[get_RLum()] 'curveType' has to be of type 'character'!")

            }

            ##RLum.type
            if (is.null(RLum.type)) {
              RLum.type <- c("RLum.Data.Curve","RLum.Data.Spectrum")

            }else if (!is(RLum.type, "character")) {
              stop("[get_RLum()] 'RLum.type' has to be of type 'character'!")

            }

            ##get.index
            if (is.null(get.index)) {
              get.index <- FALSE

            }else if (!is(get.index, "logical")) {
              stop("[get_RLum()] 'get.index' has to be of type 'logical'!")

            }

            ##get originator
            if(.hasSlot(object, "originator")){
              originator <- object@originator

            }else{
              originator <- NA_character_

            }


            ##-----------------------------------------------------------------##

            ##a pre-selection is necessary to support negative index selection
            object@records <- object@records[record.id]
            record.id <- 1:length(object@records)


            ##select curves according to the chosen parameter
            if(length(record.id)>1){

              temp <- sapply(record.id, function(x){

                if(is(object@records[[x]])[1]%in%RLum.type == TRUE){

                  ##as input a vector is allowed
                  temp <- sapply(1:length(recordType), function(k){


                    ##translate input to regular expression
                    recordType[k] <- glob2rx(recordType[k])
                    recordType[k] <- substr(recordType[k],
                                            start = 2,
                                            stop = nchar(recordType[k])-1)

                    if(grepl(recordType[k],object@records[[x]]@recordType) == TRUE &
                       object@records[[x]]@curveType%in%curveType){

                      if(!get.index){

                          object@records[[x]]

                      }else{
                        x
                      }

                    }

                  })

                  ##remove empty entries and select just one to unlist
                  temp <- temp[!sapply(temp, is.null)]

                  ##if list has length 0 skip entry
                  if(length(temp) != 0){temp[[1]]}else{temp <- NULL}

                }

              })


              ##remove empty list element
              temp <- temp[!sapply(temp, is.null)]

              ##check if the produced object is empty and show warning message
              if(length(temp) == 0){

                warning("[get_RLum()] This request produced an empty list of records!")

              }

              ##remove list for get.index
              if (get.index) {
                return(unlist(temp))

              }else{
                if (!drop) {
                  temp <- set_RLum(class = "RLum.Analysis",
                                   originator = originator,
                                   records = temp,
                                   protocol = object@protocol)
                  return(temp)

                }else{
                  if (length(temp) == 1 & recursive == TRUE) {
                    return(temp[[1]])

                  }else{
                    return(temp)

                  }

                }

              }

            }else{

              if(get.index == FALSE){


                if(drop == FALSE){

                  ##needed to keep the argument drop == TRUE
                  temp <- set_RLum(class = "RLum.Analysis",
                                   originator = originator,
                                   records = list(object@records[[record.id]]),
                                   protocol = object@protocol)
                  return(temp)

                }else{

                  return(object@records[[record.id]])

                }


              }else{

                return(record.id)

              }
            }


          })


####################################################################################################
###structure_RLum()
####################################################################################################
#' @describeIn RLum.Analysis
#' Method to show the structure of an \code{\linkS4class{RLum.Analysis}} object.
#'
#' @return
#'
#' \bold{\code{structure_RLum}}:\cr
#'
#' Returns \code{\linkS4class{data.frame}} showing the structure.
#'
#' @export
setMethod("structure_RLum",
          signature= "RLum.Analysis",
          definition = function(object) {

            ##check if the object containing other elements than allowed
            if(length(grep(FALSE, sapply(object@records, is, class="RLum.Data.Curve")))!=0){

              stop("[structure_RLum()]  Only 'RLum.Data.Curve' objects are allowed!" )

            }

            ##get length object
            temp.object.length <- length(object@records)

            ##ID
            temp.id <- 1:temp.object.length

            ##OBJECT TYPE
            temp.recordType <- c(NA)
            length(temp.recordType) <- temp.object.length
            temp.recordType <- sapply(1:temp.object.length,
                                      function(x){object@records[[x]]@recordType})

            ##PROTOCOL STEP
            temp.protocol.step <- c(NA)
            length(temp.protocol.step) <- temp.object.length

            ##n.channels
            temp.n.channels <- sapply(1:temp.object.length,
                                      function(x){length(object@records[[x]]@data[,1])})

            ##X.MIN
            temp.x.min <- sapply(1:temp.object.length,
                                 function(x){min(object@records[[x]]@data[,1])})

            ##X.MAX
            temp.x.max <- sapply(1:temp.object.length,
                                 function(x){max(object@records[[x]]@data[,1])})

            ##y.MIN
            temp.y.min <- sapply(1:temp.object.length,
                                 function(x){min(object@records[[x]]@data[,2])})

            ##X.MAX
            temp.y.max <- sapply(1:temp.object.length,
                                 function(x){max(object@records[[x]]@data[,2])})

            ##info elements as character value
            temp.info.elements <- unlist(sapply(1:temp.object.length, function(x){

              if(length(object@records[[x]]@info)!=0){
                do.call(paste, as.list(names(object@records[[x]]@info)))
              }else{NA}

            }))

            ##combine output to a data.frame
            return(
              data.frame(
                id = temp.id, recordType = temp.recordType,
                protocol.step = temp.protocol.step,
                n.channels = temp.n.channels,
                x.min = temp.x.min,
                x.max = temp.x.max,
                y.min = temp.y.min,
                y.max = temp.y.max,
                info.elements = temp.info.elements
              )
            )

          })


####################################################################################################
###length_RLum()
####################################################################################################
#' @describeIn RLum.Analysis
#' Returns the length of the object, i.e., number of stored records.
#'
#' @return
#'
#' \bold{\code{length_RLum}}\cr
#'
#' Returns the number records in this object.
#'
#' @export
setMethod("length_RLum",
          "RLum.Analysis",
          function(object){
            length(object@records)

          })

####################################################################################################
###names_RLum()
####################################################################################################
#' @describeIn RLum.Analysis
#' Returns the names of the \code{\linkS4class{RLum.Data}} objects objects (same as shown with the show method)
#'
#' @return
#'
#' \bold{\code{names_RLum}}\cr
#'
#' Returns the names of the record types (recordType) in this object.
#'
#' @export
setMethod("names_RLum",
          "RLum.Analysis",
          function(object){
            sapply(1:length(object@records), function(x){
              object@records[[x]]@recordType})

          })

