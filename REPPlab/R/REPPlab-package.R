#' Reliability Data from an Industrial Context
#' 
#' The data set contains 55 variables measured during the production
#' process of 520 units. The problem for this data is to identify the
#' produced items with a fault not detected by standard tests.
#' 
#' 
#' @name ReliabilityData
#' @docType data
#' @format A data frame with 520 observations on 55 variables.
#' @source The data was anonymized to keep confidentiality.
#' @keywords datasets
#' @examples
#' 
#' data(ReliabilityData)
#' summary(ReliabilityData) 
#' 
NULL


#' R Interface to the Java Program 'EPP-lab' v1.0
#' 
#' An interface that gives access to the Java program 'EPP-lab' which implements
#' several biologically inspired optimisation algorithms and several indices
#' for Exploratory Projection Pursuit (PP). The objective of optimizing PP
#' indices and projecting the data on the associated one-dimensional directions
#' is to detect hidden structures such as clusters or outliers in (possibly
#' high dimensional) data sets.
#' The 'EPP-lab' sources and jar-files are available under 
#' \url{https://github.com/fischuu/EPP-lab.git}. For a detailed description of
#'  'EPP-lab', see Larabi (2011).
#' 
#' \tabular{ll}{ Package: \tab REPPlab\cr Type: \tab Package\cr Version: \tab
#' 0.9.3\cr Date: \tab 2015-11-30\cr License: \tab GPL\cr LazyLoad: \tab yes\cr
#' }
#' 
#' @name REPPlab-package
#' @aliases R/EPP-lab-package
#' @docType package
#' @author Daniel Fischer, Alain Berro, Klaus Nordhausen, Anne Ruiz-Gazen
#' 
#' Maintainer: Daniel Fischer <daniel.fischer@luke.fi>
#'
#' @import rJava lattice stats LDRTools graphics
#* @importFrom utils stack
#'
#' @references \cite{Larabi Marie-Sainte, S., (2011), Biologically inspired
#' algorithms for exploratory projection pursuit, PhD thesis, University of
#' Toulouse.}
#' 
#' \cite{Ruiz-Gazen, A., Larabi Marie-Sainte, S. and Berro, A. (2010),
#' Detecting multivariate outliers using projection pursuit with particle swarm
#' optimization, \emph{COMPSTAT2010}, pp. 89-98.}
#' 
#' \cite{Berro, A., Larabi Marie-Sainte, S. and Ruiz-Gazen, A. (2010). Genetic
#' algorithms and particle swarm optimization for exploratory projection
#' pursuit. Annals of Mathematics and Artifcial Intelligence, 60, 153-178.}
#' 
#' \cite{Larabi Marie-Sainte, S., Berro, A. and Ruiz-Gazen, A. (2010). An
#' effcient optimization method for revealing local optima of projection
#' pursuit indices. \emph{Swarm Intelligence}, pp. 60-71.}
#' @keywords multivariate
NULL



