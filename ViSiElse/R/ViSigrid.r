#' Class \code{ViSigrid} defines the structure of the process to be plotted.
#' @title Class \code{ViSigrid}
#' @name ViSigrid-class
#' @rdname ViSigrid-class
#' @slot MATp A \code{"dgCMatrix"}. It stores the grid for all punctuals actions in the book. 
#' @slot MATpsup A \code{"dgCMatrix"}. It stores the grid for all punctuals actions 
#' in the book corresponding to the supplementary times.
#' @slot idsup A \code{"vector"} It stores individuals id having supplementary times. 
#' @slot colvect A \code{"matrix"} Matrix with colors to use. 
#' @slot L A \code{"data.frame"} It stores the data corresponding to long actions having a showorder. 
#' @slot idsort A \code{"matrix"} For all long actions, it stores the order of individuals in which each actions will be plot. 
#' @slot BZL A \code{"dgCMatrix"} It stores black zones for long actions, calculated for each individuals. 
#' @slot Lsup A \code{"data.frame"} for the long actions having  a showorder and suplementary times defined, it 
#' stores the data corresponding to those actions. 
#' @slot book A \code{"ViSibook"} it stores the structure of the grid for the plot. 
#' @slot group A \code{"factor"} it stores the group for the each individuals. 
#' @slot vect_tps A \code{"vector"} it stores the times vector mapping the grid. 
#' @slot informers A \code{"matrix"} It stores the indicators (mean, median or NULL) by actions. 
#' @slot testsP A \code{"vector"} Results of tests p.value<threshold.test. 
#' @slot parameters A \code{"list"}. It stores the parameters put in the \code{\link{buildViSiGrid}} function. 
#' @seealso \code{\link{buildViSiGrid}}, \code{\link{plot,ViSigrid-method}}, \code{\linkS4class{ViSibook}}
#' @exportClass ViSigrid
# Method Visigrid
ViSigrid = setClass("ViSigrid",
                   slots = c(
                     MATp = "dgCMatrix",
                     MATpsup = "dgCMatrix",
                     idsup = "vector" ,
                     colvect = "matrix",
                     L = "data.frame",
                     idsort = "matrix",
                     BZL = "dgCMatrix",
                     Lsup = "data.frame",
                     book = "ViSibook",
                     group = "factor",
                     vect_tps = "vector",
                     informers = "matrix",
                     testsP = "vector",
                     parameters = "list")
)
