#' Class \code{"TLum.BIN.File"}
#'
#' Class for luminescence curves data.
#'
#'
#' @name TLum.BIN.File-class
#' @rdname TLum.BIN.File-class
#'
#' @slot METADATA
#'  \link{data.frame}: Object containing the meta information for each curve.
#' @slot DATA
#'  \link{list}: Object containing numeric vector with count data.
#' @slot ERROR
#'   \link{list}: Object containing numeric vector with count data absolute uncertainty.
#' @slot .RESERVED
#'   \link{list}: Object containing list of undocumented raw values for internal use only.
#' @keywords classes
#'
#' @note The code and the structure of this class is based on the \linkS4class{Risoe.BINfileData} class from the \link{Luminescence} package.
#'
#' @author David Strebler
#'
#'
#' @exportClass TLum.BIN.File
#'

setClass(Class = "TLum.BIN.File",
         contains = "TLum",
         slots= c(METADATA = "data.frame",
                  DATA = "list",
                  ERROR = "list",
                  .RESERVED = "list"),
         prototype = list(METADATA = data.frame(),
                          DATA=list(),
                          ERROR=list(),
                          .RESERVED = list())
         )

# Show

#' @rdname TLum.BIN.File-class
#' @aliases show,TLum.BIN.File-method

setMethod(f = "show",
          signature = "TLum.BIN.File",
          definition = function(object){

            version <- paste(unique(object@METADATA[,"VERSION"]), collapse = ", ")
            systemID <- paste(unique(object@METADATA[,"SYSTEMID"]), collapse = ", ")
            filename <- as.character(object@METADATA[1,"FNAME"])
            nRecords <- length(object@DATA)
            records.type <- table(object@METADATA[,"LTYPE"])
            user <- paste(unique(as.character(object@METADATA[,"USER"])), collapse = ", ")
            date <- paste(unique(as.character(object@METADATA[,"DATE"])), collapse = ", ")
            run.range <- range(object@METADATA[,"RUN"])
            set.range <- range(object@METADATA[,"SET"])
            pos.range <- range(object@METADATA[,"POSITION"])

            records.type.count <- sapply(1:length(records.type),
                                         function(x){paste(
                                           names(records.type)[x],"\t(n = ",records.type[x],")",sep="")
                                         })

            records.type.count <- paste(records.type.count,
                                        collapse="\n\t                      ")

            ##print
            cat("\n[TLum.BIN.File object]")
            cat("\n\n\tBIN/BINX version     ", version)
            if(version>=6){
              cat("\n\tFile name:           ", filename)
            }
            cat("\n\t Object date:         ", date)
            cat("\n\t User:                ", user)
            cat("\n\t System ID:           ", ifelse(systemID == 0,"0 (unknown)", systemID))
            cat("\n\t Overall records:     ", nRecords)
            cat("\n\t Records type:        ", records.type.count)
            cat("\n\t Position range:      ",pos.range[1],":",pos.range[2])
            cat("\n\t Run range:           ",run.range[1],":",run.range[2])
            cat("\n\t Set range:           ",set.range[1],":",set.range[2])

          })

#set

#' @name TLum.BIN.File-class
#' @rdname TLum.BIN.File-class
#'
#' @param METADATA
#'  \link{data.frame}: Object containing the meta information for each curve.
#' @param DATA
#'  \link{list}: Object containing numeric vector with count data.
#' @param ERROR
#'   \link{list}: Object containing numeric vector with count data absolute uncertainty.
#' @param .RESERVED
#'   \link{list}: Object containing list of undocumented raw values for internal use only.
#'
#' @exportMethod set_TLum.BIN.File

setGeneric(name = "set_TLum.BIN.File",
           def = function(METADATA, DATA, ERROR, .RESERVED) {standardGeneric("set_TLum.BIN.File")}
           )

#' @rdname TLum.BIN.File-class
#' @aliases set_TLum.BIN.File set_TLum.BIN.File,TLum.BIN.File-method

setMethod(f = "set_TLum.BIN.File",
          signature = c(METADATA = "data.frame", DATA = "list", ERROR = "list", .RESERVED = "ANY"),

          definition= function(METADATA, DATA, ERROR, .RESERVED){

            if(missing(.RESERVED)){
              .RESERVED <- list()
            }

            new("TLum.BIN.File",
                METADATA = METADATA,
                DATA = DATA,
                ERROR = ERROR,
                .RESERVED = .RESERVED
            )
          })

#get

#' @name TLum.BIN.File-class
#' @rdname TLum.BIN.File-class
#'
#' @param object
#'  \linkS4class{TLum.BIN.File}: an object of class 'TLum.BIN.File'.
#'
#' @exportMethod get_TLum.BIN.File

setGeneric("get_TLum.BIN.File",
           function(object) {standardGeneric("get_TLum.BIN.File")})

#' @rdname TLum.BIN.File-class
#' @aliases get_TLum.BIN.File get_TLum.BIN.File,TLum.BIN.File-method

setMethod("get_TLum.BIN.File",
          signature=signature(object = "TLum.BIN.File"),
          definition = function(object) {

            cat("[get_TLum.FileData()]: No direct access is provided for this object type. Use the function 'TLum.BIN.File2TLum.Analysis' for object coercing.")

          })##end setMethod
