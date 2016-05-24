#' Class \code{"TLum.Analysis"}
#'
#' Object class containing analysis data for protocol analysis.
#'
#' @name TLum.Analysis-class
#' @rdname TLum.Analysis-class
#'
#' @slot protocol
#'  \link{character}: Protocol used for the analysis.
#' @slot records
#'  \link{list}: \linkS4class{TLum.Data.Curve} included in the analysis.
#'
#' @note The code and the structure of this class is based on the \linkS4class{RLum.Analysis} class from the \link{Luminescence} package.
#'
#' @keywords classes
#'
#' @author David Strebler
#'
#' @exportClass TLum.Analysis


# class definition
setClass(Class= "TLum.Analysis",
         contains = "TLum",
         slots= c(records = "list",
                  protocol = "character"),
         prototype = list (records = list(),
                           protocol = "UNKNOWN")
         )


# show method for object -------------------------------------------------------

#' @rdname TLum.Analysis-class
#' @aliases show,TLum.Analysis-method

setMethod("show",
          signature= "TLum.Analysis",
          definition=function(object){

            protocol <- object@protocol
            nRecords <- length(object@records)

            ##print
            cat("\n [TLum.Analysis]")
            cat("\n\t protocol:", protocol)
            cat("\n\t number of records:", nRecords)

            #skip this part if nothing is included in the object
            if(nRecords > 0){

              ##get object class types

              class.type <- vector()
              recordType <- vector()
              for (i in 1:nRecords){
                class.type[i] <- is(object@records[[i]])[1]
                recordType[i] <- as.character(object@records[[i]]@recordType)
              }

              table.class <- table(class.type)

              for(i in 1:length(table.class)){
                cat("\n\t .. :",names(table.class)[i],":",table.class[i])

                temp <- NULL
                k <- 1

                for(j in 1:nRecords){
                  if(names(table.class)[i] == class.type[j]){
                    temp <- paste(temp, recordType[i])
                    k <- k+1

                    if(k>10){
                      cat("\n\t .. :", temp)
                      temp <- NULL
                      k <- 1
                    }
                  }
                }

                if(!is.null(temp)){
                  cat("\n\t .. :", temp)
                  temp <- NULL
                }
              }

            }else{

              cat("\n\t >> This is an empty object and cannot be used for further analysis! <<")

            }
          }
)##end show method

# constructor (set) method for object class ------------------------------------------

#' @name TLum.Analysis-class
#' @rdname TLum.Analysis-class
#'
#' @param records
#'  \link{list}: list of \linkS4class{TLum.Data.Curve} objects
#' @param protocol
#'  \link{character}: protocol type for analysis object.
#'
#' @exportMethod set_TLum.Analysis

setGeneric("set_TLum.Analysis",
           function(records, protocol) {standardGeneric("set_TLum.Analysis")})


#' @rdname TLum.Analysis-class
#' @aliases set_TLum.Analysis set_TLum.Analysis,TLum.Analysis-method

setMethod(f = "set_TLum.Analysis",
          signature = c(records = "list",
                        protocol= "ANY"),

          definition = function(records, protocol){
            if(missing(protocol)){
              protocol <- "UNKNOWN"

            }else if (is(protocol, "character") == FALSE){
              stop("[set_TLum.Analysis] Error: 'protocol' has to be of type 'character'!")
            }

            new("TLum.Analysis",
                protocol = protocol,
                records = records
            )
          })

# constructor (get) method for object class ------------------------------------------

#' @name TLum.Analysis-class
#' @rdname TLum.Analysis-class
#'
#' @param object
#'  \linkS4class{TLum.Analysis}: an object of class TLum.Analysis.
#' @param record.id
#'  \link{numeric}: IDs of specific records.
#' @param recordType
#'  \link{character}: record type.
#' @param curveType
#'  \link{character}: curve type.
#' @param TLum.type
#'  \link{character}: TLum object type.
#' @param get.index
#'  \link{logical}: return a numeric vector with the index of each element in the TLum.Analysis object.
#' @param keep.object
#'  \link{logical}: return a TLum.Analysis object.
#'
#' @exportMethod get_TLum.Analysis

setGeneric("get_TLum.Analysis",
           function(object, record.id, recordType, curveType, TLum.type, get.index, keep.object = FALSE) {
             standardGeneric("get_TLum.Analysis")})

#' @rdname TLum.Analysis-class
#' @aliases get_TLum.Analysis get_TLum.Analysis,TLum.Analysis-method

setMethod("get_TLum.Analysis",
          signature = c(object = "TLum.Analysis",
                        record.id = "ANY",
                        recordType = "ANY",
                        curveType = "ANY",
                        TLum.type = "ANY",
                        get.index = "ANY",
                        keep.object = "ANY"),

          function(object, record.id, recordType, curveType, TLum.type, get.index, keep.object = FALSE){

            ##record.id
            if(missing(record.id)){

              record.id <- c(1:length(object@records))

            }else if (is(record.id, "numeric") == FALSE){

              stop("[get_TLum.Analysis()] 'record.id' has to be of type 'numeric'!")

            }

            ##check if record.id exists
            if(FALSE%in%(abs(record.id)%in%(1:length(object@records)))){

              stop("[get_TLum.Analysis()] At least one 'record.id' is invalid!")

            }

            ##recordType
            if(missing(recordType)){

              recordType <- unique(
                unlist(
                  lapply(1:length(object@records),
                         function(x){object@records[[x]]@recordType})))

            }else{

              if (is(recordType, "character") == FALSE){

                stop("[get_TLum.Analysis()] Error: 'recordType' has to be of type 'character'!")

              }

            }

            ##curveType
            if(missing(curveType) == TRUE){

              curveType <- unique(
                unlist(
                  lapply(1:length(object@records),
                         function(x){object@records[[x]]@curveType})))

            }else if (is(curveType, "character") == FALSE){

              stop("[get_TLum.Analysis()] Error: 'curveType' has to be of type 'character'!")

            }

            ##TLum.type
            if(missing(TLum.type) == TRUE){

              TLum.type <- c("TLum.Data.Curve","TLum.Data.Spectrum")

            }else if (is(TLum.type, "character") == FALSE){

              stop("[get_TLum.Analysis()] Error: 'TLum.type' has to be of type 'character'!")

            }

            ##get.index
            if(missing(get.index) == TRUE){

              get.index <- FALSE

            }else if (is(get.index, "logical") == FALSE){

              stop("[get_TLum.Analysis()] Error: 'get.index' has to be of type 'logical'!")

            }




            ##-----------------------------------------------------------------##

            ##a pre-selection is necessary to support negative index selection
            object@records <- object@records[record.id]
            record.id <- 1:length(object@records)


            ##select curves according to the chosen parameter
            if(length(record.id)>1){

              temp <- sapply(record.id, function(x){

                if(is(object@records[[x]])[1]%in%TLum.type == TRUE){

                  ##as input a vector is allowed
                  temp <- sapply(1:length(recordType), function(k){


                    ##translate input to regular expression
                    recordType[k] <- glob2rx(recordType[k])
                    recordType[k] <- substr(recordType[k],
                                            start = 2,
                                            stop = nchar(recordType[k])-1)

                    if(grepl(recordType[k],object@records[[x]]@recordType) == TRUE &
                       object@records[[x]]@curveType%in%curveType){

                      if(get.index == FALSE){

                        object@records[[x]]

                      }else{x}

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

              ##check if produced object is empty and show warning message
              if(length(temp) == 0){

                warning("This request has produced an empty 'TLum.Analysis' object!")

              }

              ##remove list for get.index
              if(get.index == TRUE){

                return(unlist(temp))

              }else{

                if(keep.object == TRUE){

                  temp <- set_TLum.Analysis(records = temp, protocol = object@protocol)
                  return(temp)

                }else{

                  if(length(temp) == 1){

                    return(temp[[1]])

                  }else{

                    return(temp)

                  }

                }

              }

            }else{

              if(get.index == FALSE){


                if(keep.object == TRUE){

                  ##needed to keep the argument keep.object == TRUE
                  temp <- set_TLum.Analysis(records = list(object@records[[record.id]]),
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
