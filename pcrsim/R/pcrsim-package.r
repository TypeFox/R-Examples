###############################################################################
#' Forensic DNA process simulator.
#'
#' PCRsim is a package to simulate the forensic DNA process. The function \code{pcrsim}
#' opens up a graphical user interface which allow the user to enter parameters
#' required for the simulation. Once calibrated the program can potentially
#' be used to: reduce the laboratory work needed to validate new STR kits,
#' create samples for educational purposes, help develop methods for
#' interpretation of DNA evidence, etc.
#' 
#' This is a first version which is still experimental and under development.
#' 
#' Areas in need of more research are better calibration and more correct
#' scaling to peak heights over a range of input amounts. The current
#' implementation is built to mimic the biological processes as closely as
#' possible and are not suitable for simulation of large number of samples
#' due to performance.
#' 
#' @title Simulation of the Forensic DNA process
#' @docType package
#' @name pcrsim-package
#' @author Oskar Hansson \email{oskar.hansson@@fhi.no}
#' @keywords package
#' @references  Gill, Peter, James Curran, and Keith Elliot.
#' \\u0022 A Graphical Simulation Model of the Entire DNA Process Associated with
#'  the Analysis of Short Tandem Repeat Loci\\u0022
#'  Nucleic Acids Research 33, no. 2 (2005): 632-643. doi:10.1093/nar/gki205.
#'   
NULL