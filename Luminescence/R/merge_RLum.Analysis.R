#' Merge function for RLum.Analysis S4 class objects
#'
#' Function allows merging of RLum.Analysis objects and adding of allowed
#' objects to an RLum.Analysis.
#'
#' This function simply allowing to merge \code{\linkS4class{RLum.Analysis}}
#' objects.  Additionally other \code{\linkS4class{RLum}} objects can be added
#' to an existing \code{\linkS4class{RLum.Analysis}} object. Supported objects
#' to be added are: \code{\linkS4class{RLum.Data.Curve}},
#' \code{\linkS4class{RLum.Data.Spectrum}} and
#' \code{\linkS4class{RLum.Data.Image}}.\cr
#'
#' The order in the new \code{\linkS4class{RLum.Analysis}} object is the object
#' order provided with the input list.
#'
#' @param objects \code{\link{list}} of \code{\linkS4class{RLum.Analysis}}
#' (\bold{required}): list of S4 objects of class \code{RLum.Analysis}.
#' Furthermore other objects of class \code{\linkS4class{RLum}} can be added,
#' see details.
#' @return Return an \code{\linkS4class{RLum.Analysis}} object.
#' @note The information for the slot 'protocol' is taken from the first
#' \code{\linkS4class{RLum.Analysis}} object in the input list. Therefore at
#' least one object of type \code{\linkS4class{RLum.Analysis}} has to be
#' provided.
#' @section Function version: 0.1
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#' @seealso \code{\link{merge_RLum}}, \code{\linkS4class{RLum.Analysis}},
#' \code{\linkS4class{RLum.Data.Curve}},
#' \code{\linkS4class{RLum.Data.Spectrum}},
#' \code{\linkS4class{RLum.Data.Image}}, \code{\linkS4class{RLum}}
#' @references -
#' @keywords utilities
#' @examples
#'
#'
#' ##merge different RLum objects from the example data
#' data(ExampleData.RLum.Analysis, envir = environment())
#' data(ExampleData.BINfileData, envir = environment())
#'
#' object <- Risoe.BINfileData2RLum.Analysis(CWOSL.SAR.Data, pos=1)
#' curve <- get_RLum(object)[[2]]
#'
#' temp.merged <- merge_RLum.Analysis(list(curve, IRSAR.RF.Data, IRSAR.RF.Data))
#'
#' @export
merge_RLum.Analysis<- function(
  objects
){

  # Ingegrity checks ----------------------------------------------------------------------------

  ##check if object is of class RLum
  temp.class.test <- sapply(1:length(objects), function(x){

    if(is(objects[[x]], "RLum") == FALSE){

      temp.text <- paste("[merge_RLum.Analysis()]: At least element", x, "is not of class 'RLum' or a derivative class!")
      stop(temp.text)
    }




    ##provide class of objects
    is(objects[[x]])[1]

  })

  ##check if at least one object of RLum.Analysis is provided
  if(!"RLum.Analysis"%in%temp.class.test){

    stop("[merge_RLum.Analysis()] At least one input object in the list has to be of class
           'RLum.Analysis'!")

  }



  # Merge objects -------------------------------------------------------------------------------

  ##(0) get recent environment to later set variable temp.meta.data.first
  temp.environment  <- environment()
  temp.meta.data.first <- NA; rm(temp.meta.data.first) #to avoid problems with the R check routine

  ##(1) collect all elements a a list
  temp.element.list <- unlist(sapply(1:length(objects), function(x){

    ##Depending on the element the right functions is used
    if(is(objects[[x]])[1] == "RLum.Analysis"){

      ##grep export meta data from the first RLum.Analysis objects an write
      if(!exists("temp.meta.data.first")){

        assign("temp.meta.data.first", objects[[x]]@protocol, envir = temp.environment)

      }

      ##return to list
      get_RLum(objects[[x]])

    }else if((is(objects[[x]])[1] == "RLum.Data.Curve") |
               (is(objects[[x]])[1] == "RLum.Data.Image") |
               (is(objects[[x]])[1] == "RLum.Data.Spectrum")){

      ##return to list
      objects[[x]]

    }else{

      stop("[merge_RLum.Anlysis()] What ever you provided, this 'RLum' object is not supported here!")

    }


  }))


  # Build new RLum.Analysis object --------------------------------------------------------------

  temp.new.RLum.Analysis <- set_RLum(
    class = "RLum.Analysis",
    originator = "merge_RLum.Analysis",
    records = temp.element.list,
    protocol = temp.meta.data.first)


  # Return object -------------------------------------------------------------------------------

  return( temp.new.RLum.Analysis)

}
