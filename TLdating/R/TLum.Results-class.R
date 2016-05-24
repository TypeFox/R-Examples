#' Class \code{"TLum.Results"}
#'
#' Object class contains results data from functions.
#'
#'
#' @name TLum.Results-class
#' @rdname TLum.Results-class
#'
#' @slot originator
#'  \link{character}: contains the name of the calling function (the function that produces this object).
#' @slot data
#'  \link{list}: a list containing the data to be stored in the object.
#'
#' @note The code and the structure of this class is based on the \linkS4class{RLum.Results} class from the \link{Luminescence} package.
#'
#' @keywords classes
#'
#' @author David Strebler
#'
#' @exportClass TLum.Results


##class definition
setClass(Class="TLum.Results",
         contains = "TLum",
         slots = c(originator = "character",
                   data = "list"),
         prototype = list(originator = character(),
                          data = list())
         )


# show method for object ------------------------------------------------------

#' @rdname TLum.Results-class
#' @aliases show,TLum.Results-method

setMethod("show",
          signature(object = "TLum.Results"),
          function(object){


            ##data elements
            temp.names <- names(object@data)
            temp.type <- sapply(1:length(object@data),
                                function(x){

                                  paste("\t .. $", temp.names[x],
                                        " : ",
                                        is(object@data[[x]])[1],
                                        sep = "")


                                })
            temp.type <- paste(temp.type, collapse="\n")

            ##print information
            cat("\n [TLum.Results]")
            cat("\n\t originator: ", object@originator,"()", sep="")
            cat("\n\t data:", length(object@data))
            cat("\n", temp.type)


          }
)



# constructor (set) method for object class -------------------------------

#' @name TLum.Results-class
#' @rdname TLum.Results-class
#'
#' @param originator
#'  \link{character}: : contains the name of the calling function.
#' @param data
#'  \link{list}:  the data to be stored in the object.
#'
#' @exportMethod set_TLum.Results

setGeneric("set_TLum.Results",
           function(originator, data) {standardGeneric("set_TLum.Results")})

#' @rdname TLum.Results-class
#' @aliases set_TLum.Results set_TLum.Results,TLum.Results-method

setMethod(f = "set_TLum.Results",
          signature = c(originator = "ANY",
                        data = "list"),

          definition = function(originator, data){

            if(missing(originator)){
              originator <- "Unknown"

            }else if(!is.character(originator)){
              stop("[set_TLum.Results] Error: 'originator' is not of type 'character'.")
            }

            if(missing(data)){
              stop("[set_TLum.Results] Error: 'data' is missing.")
            }

            new("TLum.Results",
                originator = originator,
                data = data)
          })


# GetMethods --------------------------------------------------------------

#' @name TLum.Results-class
#' @rdname TLum.Results-class
#'
#' @param object
#'  \linkS4class{TLum.Results}:  object to be evaluated.
#' @param ref
#'  \link{character}: name of the 'data' slot to be returned.
#'
#' @exportMethod get_TLum.Results

setGeneric("get_TLum.Results",
           function(object, ref) {standardGeneric("get_TLum.Results")})

#' @rdname TLum.Results-class
#' @aliases get_TLum.Results get_TLum.Results,TLum.Results-method

setMethod("get_TLum.Results",
          signature=signature(object = "TLum.Results",
                              ref = "ANY"),
          definition = function(object, ref) {

            if(missing(ref)){
              res <- object@data

            }else if(!is.character(ref)){
                stop("[get_TLum.Results] Error: 'ref' has to be a character!")

            }else if(ref %in% names(object@data)){
              res <- object@data[[ref]]

            }else if(ref == "data"){
              res <- object@data

            }else if(ref == "originator"){
              res <- object@originator

            }else{
              stop("[get_TLum.Results] Error: Input 'ref' is unknown.")
            }

            return(res)
          })
