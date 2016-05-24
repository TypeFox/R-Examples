#' @title Data to support the FSA package.
#' 
#' @description This package contains data to support the FSA package.
#' 
#' @details This package contains additional data files that can be used for common fisheries stock assessment methods described in the \pkg{FSA} package and on the \href{http://derekogle.com/fishR/}{fishR website}.
#' 
#' The help files for these datasets are embedded with topics that can be searched to find data files that can be analyzed with those topics.  For example, use the following commands to find data files for the corresponding topics.
#' 
#' \tabular{ll}{
#'  \code{help.search("Length Expansion",package=c("FSAdata","FSA"))} \tab Expand subsampled lengths.\cr
#'  \code{help.search("Length Conversion",package=c("FSAdata","FSA"))} \tab Convert between length types.\cr
#'  \code{help.search("Age Comparison",package=c("FSAdata","FSA"))} \tab Ageing (error, precision, or comparison).\cr
#'  \code{help.search("Age-Length Key",package=c("FSAdata","FSA"))} \tab Age-Length Key data.\cr
#'  \code{help.search("Weight-Length",package=c("FSAdata","FSA"))} \tab Weight-length model data.\cr
#'  \code{help.search("Length Frequency",package=c("FSAdata","FSA"))} \tab Length frequency data.\cr
#'  \code{help.search("Size Structure",package=c("FSAdata","FSA"))} \tab Size structure data.\cr
#'  \code{help.search("Abundance",package=c("FSAdata","FSA"))} \tab Data for abundance estimates.\cr
#'  \code{help.search("Capture-Recapture",package=c("FSAdata","FSA"))} \tab Mark-recapture data.\cr
#'  \code{help.search("Mark-Recapture",package=c("FSAdata","FSA"))} \tab Mark-recapture data.\cr
#'  \code{help.search("Capture History",package=c("FSAdata","FSA"))} \tab Capture history mark-recapture (compare to summarized data) data.\cr
#'  \code{help.search("Petersen",package=c("FSAdata","FSA"))} \tab Petersen mark-recapture (closed population, single sample).\cr
#'  \code{help.search("Schnabel",package=c("FSAdata","FSA"))} \tab Schnabel mark-recapture (closed population, multiple samples).\cr
#'  \code{help.search("Jolly-Seber",package=c("FSAdata","FSA"))} \tab Jolly-Seber mark-recapture (open population, multiple samples).\cr
#'  \code{help.search("Depletion",package=c("FSAdata","FSA"))} \tab Depletion (Leslie, DeLury) methods for estimating abundance.\cr
#'  \code{help.search("Removal",package=c("FSAdata","FSA"))} \tab Removal (K-pass) methods for estimating abundance.\cr 
#'  \code{help.search("Mortality",package=c("FSAdata","FSA"))} \tab Data for mortality estimation.\cr
#'  \code{help.search("Catch curve",package=c("FSAdata","FSA"))} \tab Catch curve.\cr

#'  \code{help.search("Growth",package=c("FSAdata","FSA"))} \tab Growth model data.\cr
#'  \code{help.search("Recruitment",package=c("FSAdata","FSA"))} \tab Stock-recruitment and recruitment time-series data.\cr
#'  \code{help.search("Maturity",package=c("FSAdata","FSA"))} \tab Maturity data.\cr
#' }
#'  
#' Additional fisheries-related data sets are in the \code{\link[FSA]{FSA}} and \code{fishmethods} packages.
#' 
#' @docType package
#' 
#' @name FSAdata
#' 
#' @export
FSAdataTopics <- c("Length Expansion","Length Conversion",
                   "Age Comparison","Age-Length Key","Back-Calculation",
                   "Weight-Length","Length Frequency","Size Structure",
                   "Capture-Recapture","Depletion","Removal",
                   "Mortality","Growth","Recruitment","Maturity",
                   "Other")