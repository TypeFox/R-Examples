#' General plot function for RLum S4 class objects
#'
#' Function calls object specific plot functions for RLum S4 class objects.
#'
#' The function provides a generalised access point for plotting specific
#' \code{\linkS4class{RLum}} objects.\cr Depending on the input object, the
#' corresponding plot function will be selected.  Allowed arguments can be
#' found in the documentations of each plot function.  \tabular{lll}{
#' \bold{object} \tab \tab \bold{corresponding plot function} \cr
#'
#' \code{\linkS4class{RLum.Data.Curve}} \tab : \tab
#' \code{\link{plot_RLum.Data.Curve}} \cr
#' \code{\linkS4class{RLum.Data.Spectrum}} \tab : \tab
#' \code{\link{plot_RLum.Data.Spectrum}}\cr
#' \code{\linkS4class{RLum.Data.Image}} \tab : \tab
#' \code{\link{plot_RLum.Data.Image}}\cr \code{\linkS4class{RLum.Analysis}}
#' \tab : \tab \code{\link{plot_RLum.Analysis}}\cr
#' \code{\linkS4class{RLum.Results}} \tab : \tab
#' \code{\link{plot_RLum.Results}} }
#'
#' @param object \code{\linkS4class{RLum}} (\bold{required}): S4 object of
#' class \code{RLum}. Optional a \code{\link{list}} containing objects of class \code{\linkS4class{RLum}}
#' can be provided. In this case the function tries to plot every object in this list according
#' to its \code{RLum} class.
#'
#' @param \dots further arguments and graphical parameters that will be passed
#' to the specific plot functions. The only argument that is supported directly is \code{main}
#' (setting the plot title). In contrast to the normal behaviour \code{main} can be here provided as
#' \code{\link{list}} and the arguments in the list will dispatched to the plots if the \code{object}
#' is of type \code{list} as well.
#'
#' @return Returns a plot.
#'
#' @note The provided plot output depends on the input object.
#'
#' @section Function version: 0.4.2
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#'
#' @seealso \code{\link{plot_RLum.Data.Curve}},
#' \code{\linkS4class{RLum.Data.Curve}}, \code{\link{plot_RLum.Data.Spectrum}},
#' \code{\linkS4class{RLum.Data.Spectrum}}, \code{\link{plot_RLum.Data.Image}},
#' \code{\linkS4class{RLum.Data.Image}}, \code{\link{plot_RLum.Analysis}},
#' \code{\linkS4class{RLum.Analysis}}, \code{\link{plot_RLum.Results}},
#' \code{\linkS4class{RLum.Results}}
#'
#' @references #
#'
#' @keywords dplot
#'
#' @examples
#'
#'
#' #load Example data
#' data(ExampleData.CW_OSL_Curve, envir = environment())
#'
#' #transform data.frame to RLum.Data.Curve object
#' temp <- as(ExampleData.CW_OSL_Curve, "RLum.Data.Curve")
#'
#' #plot RLum object
#' plot_RLum(temp)
#'
#'
#' @export
plot_RLum<- function(
  object,
  ...
){

  # Define dispatcher function ----------------------------------------------------------

  ##check if object is of class RLum
  RLum.dispatcher <- function(object, ...) {
    if (inherits(object, "RLum")) {

      ##grep object class
      object.class <- is(object)[1]

      ##select which plot function should be used and call it
      switch (
        object.class,

        RLum.Data.Curve = plot_RLum.Data.Curve(object = object, ...),
        RLum.Data.Spectrum = plot_RLum.Data.Spectrum(object = object, ...),
        RLum.Data.Image = plot_RLum.Data.Image(object = object, ...),
        RLum.Analysis = plot_RLum.Analysis(object = object, ...),
        RLum.Results = plot_RLum.Results(object = object, ...)

      )

    }else{
      stop(paste0(
        "[plot_RLum()] Sorry, I don't know what to do for object of type '", is(object)[1], "'."
      ))

    }

  }



  # Run dispatcher ------------------------------------------------------------------------------

  ##call for the list, if not just proceed as normal
  if(is(object, "list")) {
    ##(1) get rid of objects which are not RLum objects to avoid errors
    object.cleaned <-
      object[sapply(object, inherits, what = "RLum")]

    ##(1.1) place warning message
    if (length(object) > length(object.cleaned)) {
      warning(paste0(
        length(object) - length(object.cleaned)," non 'RLum' object(s) removed from list."
      ))

    }

    ##(2) check if empty, if empty do nothing ...
    if (length(object.cleaned) != 0) {

      ## If we iterate over a list, this might be extremly useful to have different plot titles
      if("main" %in% names(list(...))){
        if(is(list(...)$main,"list")){
          main.list <- rep(list(...)$main, length = length(object.cleaned))

        }
      }

      ##set also mtext, but in a different way
      if(!"mtext" %in% names(list(...))){


        if(is(object[[1]], "RLum.Analysis")){
          mtext <- paste("Record:", 1:length(object.cleaned))

        }else{
          mtext <- NULL

        }
      }else{
        mtext <- rep(list(...)$mtext, lenth.out = length(object.cleaned))

      }


      if(exists("main.list")){
        ##dispatch objects
        for (i in 1:length(object.cleaned)) {
          RLum.dispatcher(object = object[[i]],
                          main = main.list[[i]],
                          mtext = mtext[[i]],
                          ...)
        }
      }else{
        for (i in 1:length(object.cleaned)) {


          RLum.dispatcher(object = object[[i]],
                          mtext = mtext[[i]],
                          ...)
        }

      }

    }

  }else{
    ##dispatch object
    RLum.dispatcher(object = object, ...)

  }

}

