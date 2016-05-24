#' MPI Profiling Tools
#' 
#' This package contains several libraries for profiling MPI codes, as well as
#' some R tools for parsing, analyzing, and plotting their outputs.
#' 
#' \tabular{ll}{ Package: \tab pbdPROF\cr Type: \tab Package\cr License: \tab
#' MPL\cr LazyLoad: \tab yes\cr } This package requires an MPI library
#' (OpenMPI, MPICH2, or LAM/MPI).
#' 
#' @docType package
#' @author Wei-Chen Chen, Drew Schmidt, Gaurav Sehrawat, Pragneshkumar Patel,
#' and George Ostrouchov.
#' @references Programming with Big Data in R Website: \url{http://r-pbd.org/}
#' @keywords Package
#' @examples
#' 
#' \dontrun{
#' demo(allreduce.fpmpi, "pbdPROF")
#' demo(svd.fpmpi, "pbdPROF")
#' demo(masterslavePI.fpmpi, "pbdPROF")
#' demo(allreduce.mpip, "pbdPROF")
#' demo(svd.mpip, "pbdPROF")
#' demo(masterslavePI.mpip, "pbdPROF")
#' }
#' 
#' @name pbdPROF-package
#' @aliases pbdPROF pbdPROF-package pbdPROF
NULL





#' Profilers
#' 
#' Profilers Directly Supported by pbdPROF
#' 
#' 
#' @section Details: Currently, the \pkg{fpmpi} library is fully supported; a
#' version of this library is bundled with the source of the pbdPROF package
#' (see package vignette for more details).  Additional libraries can easily be
#' linked with pbdPROF, but these are not yet fully supported.  The \pkg{mpiP}
#' and \pkg{TAU} profilers are expected to be fully supported by the conclusion
#' of Google Summer of Code (~September 30 2013).
#' @references Programming with Big Data in R Website: \url{http://r-pbd.org/}
#' 
#' \pkg{fpmpi} website:
#' \url{http://www.mcs.anl.gov/research/projects/fpmpi/WWW/}
#' 
#' \pkg{mpiP} website: \url{http://mpip.sourceforge.net/}
#' 
#' \pkg{TAU} website: \url{http://www.cs.uoregon.edu/research/tau/home.php}
#' 
#' @name profilers
NULL


