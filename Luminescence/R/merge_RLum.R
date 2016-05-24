#' General merge function for RLum S4 class objects
#'
#' Function calls object-specific merge functions for RLum S4 class objects.
#'
#' The function provides a generalised access point for merge specific
#' \code{\linkS4class{RLum}} objects.\cr Depending on the input object, the
#' corresponding merge function will be selected.  Allowed arguments can be
#' found in the documentations of each merge function. Empty list elements (\code{NULL}) are
#' automatically removed from the input \code{list}.
#'
#' \tabular{lll}{
#' \bold{object} \tab \tab \bold{corresponding merge function} \cr
#'
#' \code{\linkS4class{RLum.Results}} \tab : \tab
#' \code{merge_RLum} }
#'
#' @param objects \code{\link{list}} of \code{\linkS4class{RLum}}
#' (\bold{required}): list of S4 object of class \code{RLum}
#'
#' @param \dots further arguments that one might want to pass to the specific
#' merge function
#'
#' @return Return is the same as input objects as provided in the list.
#'
#' @note So far not for every \code{RLum} object a merging function exists.
#'
#' @section Function version: 0.1.2
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#'
#' @seealso \code{\linkS4class{RLum.Data.Curve}}, \code{\linkS4class{RLum.Data.Image}},
#' \code{\linkS4class{RLum.Data.Spectrum}}, \code{\linkS4class{RLum.Analysis}}, \code{\linkS4class{RLum.Results}}
#'
#' @references #
#'
#' @keywords utilities
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
#' ##apply the central dose model 2nd time
#' temp2 <- calc_CentralDose(ExampleData.DeValues$CA1)
#'
#' ##merge the results and store them in a new object
#' temp.merged <- get_RLum(merge_RLum(objects = list(temp1, temp2)))
#'
#'
#' @export
merge_RLum<- function(
  objects,
  ...
){

  # Integrity check ----------------------------------------------------------
    if(!is.list(objects)){
      stop("[merge_RLum()] argument 'objects' needs to be of type list!")

    }

    ##we are friendly and remove all empty list elements, this helps a lot if we place things
    ##we DO NOT provide a warning as this lower the computation speed in particular cases.
    objects <- objects[!sapply(objects, is.null)]

  ##if list is empty afterwards we do nothing
   if(length(objects) != 0) {
      ##check if objects are of class RLum
      temp.class.test <- unique(sapply(1:length(objects), function(x) {
        if (!is(objects[[x]], "RLum")) {
          temp.text <-
            paste(
              "[merge_RLum()]: At least element", x, "is not of class 'RLum' or a derivative class!"
            )
          stop(temp.text, call. = FALSE)
        }
        ##provide class of objects ... so far they should be similar
        is(objects[[x]])[1]
      }))

      ##check if objects are consitent
      if (length(temp.class.test) > 1) {
        ##This is not valid for RLum.Analysis objects
        if (!"RLum.Analysis" %in% temp.class.test) {
          stop("[merge_RLum()] So far only similar input objects in the list are supported!")
        }
      }

      ##grep object class
      objects.class <-
        ifelse("RLum.Analysis" %in% temp.class.test, "RLum.Analysis", temp.class.test)

      ##select which merge function should be used
      switch (
        objects.class,
        RLum.Data.Image = stop(
          "[merge_RLum()] Sorry, merging of 'RLum.Data.Image' objects is currently not supported!"
        ),
        RLum.Data.Spectrum = stop(
          "[merge_RLum()] Sorry, merging of 'RLum.Data.Spectrum' objects is currently not supported!"
        ),
        RLum.Data.Curve = merge_RLum.Data.Curve(objects, ...),
        RLum.Analysis = merge_RLum.Analysis(objects, ...),
        RLum.Results = merge_RLum.Results(objects, ...)
      )

    }else{

      warning("[merge_RLum()] Nothing was merged as object list was found to be empty!")
      return(NULL)

    }

}
