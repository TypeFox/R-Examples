#' Bundle of functions for modifying MAESTRA/MAESPA input files, and reading
#' output files.
#' 
#' The main functions are \code{\link{runmaespa}} and
#' \code{\link{maesparunall}}, see their help pages. Functions that read
#' parameters, and modify parameters or whole namelists are
#' \code{\link{readPAR}}, \code{\link{replacePAR}},
#' \code{\link{replaceNameList}}, and \code{\link{parseFile}} to read an entire
#' file. Functions that read output are \code{\link{readdayflux}},
#' \code{\link{readhrflux}}
#' 
#' \tabular{ll}{ Package: \tab Maeswrap\cr Type: \tab Package\cr Version: \tab
#' 1.2\cr Date: \tab 2008-12-03\cr License: \tab GPL\cr LazyLoad: \tab yes\cr }
#' 
#' @name Maeswrap-package
#' @aliases Maeswrap-package Maeswrap
#' @docType package
#' @author Remko Duursma Maintainer: Remko Duursma <remkoduursma@@gmail.com>
#' @references See Belinda Medlyn's MAESTRA homepage at:
#' \url{http://maespa.github.io}
#' @keywords package
NULL





#' Example Maeswrap definition file
#' 
#' The MAESTRA/MAESPA wrapper needs a 'definition file', where names of
#' parameters are defined, together with their locations in parameter files and
#' namelists. The 'comment' column is ignored when parsing the definition file.
#' 
#' 
#' @name maeswrapdefinitions
#' @docType data
#' @format A white space separated dataset (readable with \code{read.table}).
#' @references None.
#' @keywords datasets
NULL





#' Example Maeswrap run datafile.
#' 
#' Example of a 'run dataset', where each row corresponds to a simulation with
#' parameters set by the column names. The wrapper looks for the column names
#' in the definition file, and sets the parameters based on information in that
#' file.
#' 
#' 
#' @name runfiletest
#' @docType data
#' @format A comma-separated file.
#' @references None.
#' @keywords datasets
NULL



