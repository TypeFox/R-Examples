#' Microbead Surface Nucleic Acid Melting Curve Analysis
#' 
#' Nucleic acid Melting Curve Analysis is a powerful method to 
#' investigate the interaction of double stranded nucleic acids. The MBmca package 
#' provides data sets and lightweight utilities for nucleic acid melting curve 
#' analysis and presentation on microbead surfaces. However, the function of the 
#' package can also be used for the analysis of reactions in solution (e.g., qPCR). 
#' Methods include melting curve data pre-processing (smooth, normalize, rotate, 
#' background subtraction), data inspection (comparison of multiplex melting 
#' curves) with location parameters (mean, median), deviation parameters (standard 
#' of the melting peaks including the second derivative. The second derivative 
#' melting peaks is implemented as parameter to further characterize the melting 
#' behavior. Plot functions to illustrate data quality, smoothed curves and 
#' derivatives are available too.
#' 
#' \tabular{ll}{ Package: \tab MBmca\cr Type: \tab Package\cr Version: \tab
#' 0.0.3-4\cr Date: \tab 2014-08-17\cr License: GPL (>= 2) }
#' 
#' @name MBmca-package
#' @aliases MBmca-package MBmca
#' @docType package
#' @author Stefan Roediger <stefan_roediger@@gmx.de>
#' @references A Highly Versatile Microscope Imaging Technology Platform for
#' the Multiplex Real-Time Detection of Biomolecules and Autoimmune Antibodies.
#' S. Roediger, P. Schierack, A. Boehm, J. Nitschke, I. Berger, U. Froemmel, C.
#' Schmidt, M. Ruhland, I. Schimke, D. Roggenbuck, W. Lehmann and C.
#' Schroeder.  \emph{Advances in Biochemical Bioengineering/Biotechnology}.
#' 133:35--74, 2013. \url{http://www.ncbi.nlm.nih.gov/pubmed/22437246}
#' 
#' Surface Melting Curve Analysis with R. S. Roediger, A. Boehm and I.
#' Schimke. \emph{The R Journal}. 5(2):37--52, 2013.
#' \url{http://journal.r-project.org}
#' 
#' Nucleic acid detection based on the use of microbeads: a review. S.
#' Roediger, C. Liebsch, C. Schmidt, W. Lehmann, U. Resch-Genger, U. Schedler,
#' P. Schierack. \emph{Microchim Acta} 2014:1--18. DOI:
#' 10.1007/s00604-014-1243-4
#' @keywords package
#' @importFrom chipPCR inder
#' @importFrom robustbase lmrob
#' 
#' 
NULL









