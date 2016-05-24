#' @name miscset-package
#' @docType package
#' @useDynLib miscset
#' @aliases miscset
#' @keywords misc miscellaneous tools
#' @author Sven E. Templer
#' @title Miscellaneous R Tools
#' @description 
#' A collection of miscellaneous methods to simplify various tasks,
#' including plotting, data.frame and matrix transformations, environment
#' functions, regular expression methods, and string and logical operations, as
#' well as numerical and statistical tools. Most of the methods are simple but
#' useful wrappers of common base R functions, which extend S3 generics or
#' provide default values for important parameters.
#' @details 
#' The package vignette provides a comprehensive overview and examples for the
#' usage of all available functions in the package. View with 
#' \code{vignette("miscset")}.

NULL

#' @importFrom parallel mclapply
#' @importFrom xtable xtable print.xtable
#' @importFrom gridExtra grid.arrange
#' @importFrom tools startDynamicHelp
#' @importFrom ggplot2 ggplot geom_point theme_bw theme element_rect 
#' scale_x_discrete scale_y_discrete
#' @importFrom grDevices cairo_pdf cairo_ps dev.off hcl pdf png postscript 
#' setEPS svg
#' @importFrom graphics arrows barplot plot
#' @importFrom stats confint qnorm sd symnum
#' @importFrom utils capture.output combn help object.size packageVersion 
#' tail

.onAttach <- function(libname, pkgname) {
  if (interactive()) { packageStartupMessage(
    'miscset version ',
    as.character(packageVersion("miscset")),
    ' | help and features: vignette("miscset")')
  } }
