#' knitLatex: Latex table helpers for knitr.
#'
#' This package was inspired by the xtable package, but allows for more
#' fine_grained control, especially in regards to the longtable and supertabular
#' (which is not included in xtable) environments. This package provides four
#' functions to assist in using knitr with latex:
#'
#' \itemize{
#'   \item xTab: creates a basic latex table.
#'   \item lTab: creates a longtable evironment.
#'   \item sTab: creates a supertabular environment.
#'   \item knitr_sethooks: fixes a bug in the knit_hook chunk and provids a
#'     'com' hook which turns knitr output into latex commands
#' }
#'
#' @docType package @name knitLatex
NULL
