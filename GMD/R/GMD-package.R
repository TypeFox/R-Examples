##' Computate Generalized Minimum Distance (GMD) between
##' discrete distributions and clustering tools
##'
##' \tabular{ll}{
##' Package: \tab GMD \cr
##' Type: \tab Package \cr
##' License: \tab GPL (>= 2) \cr
##' }
##' 
##' This package contains: \cr
##' 1) modules and functions for GMD computation, with GMD algorithm implemented in C to interface with R.\cr
##' 2) related clustering and visualization tools for distributions.\cr
##' \cr
##' 
##' An overview of functions
##' \tabular{ll}{
##' Function \tab Description\cr
##' \code{ghist} \tab  Generalized Histogram Computation and Visualization\cr
##' \code{gdist} \tab  Generalized Distance Matrix Computation\cr
##' \code{css} \tab  Computing Clustering Sum-of-Squares and\cr
##'              \tab evaluating the clustering by the \emph{``elbow''} method\cr
##' \code{heatmap.3} \tab  Enhanced Heatmap Representation with Dendrogram and Partition\cr
##' \code{gmdp} \tab Computation of GMD on a pair of histograms \cr
##' \code{gmdm} \tab Computation of GMD Matrix on a set of histograms \cr
##' }
##' 
##'
##' ## To install from online repositories (e.g. CRAN)
##' install.packages(pkgs="GMD", repos="http://cran.r-project.org")
##' 
##' ## Sometimes the offical repository might not be up to date, then
##' ## you may install it from a downloaded source file; please replace
##' ## '<current-version>' with actual version numbers: Note that as
##' ## new versions are release, the '<current-version>' changes.
##' install.packages(pkgs="GMD_<current-version>.tar.gz", repos=NULL)
##' 
##' ## Load the package and get a complete list of functions, use
##' library(GMD)
##' ls("package:GMD")
##'
##' ## help documantation of the package
##' help(GMD)  # this page
##' 
##' @name GMD-package
##' @aliases GMD GMD-package
##' @docType package
##' @title The Package for Generalized Minimum Distance (GMD) Computation
##' @references Zhao et al (2011),
##' "Systematic Clustering of Transcription Start Site Landscapes",
##' \emph{PLoS ONE} \bold{6}(8): e23409.\cr
##' \url{http://dx.plos.org/10.1371/journal.pone.0023409}
##'
##' See \code{citation("GMD")} for BibTeX entries for LaTeX users.
##' @author Xiaobei Zhao and Albin Sandelin
##' 
##' Maintainer: Xiaobei Zhao \email{xiaobei (at) binf.ku.dk}
##' @concept GMD gmd distance nonparametric optimize cluster classif
##' @keywords package
##' @return \code{NULL}
##' @seealso \code{\link{gmdp}}, \code{\link{gmdm}}, \code{\link{cage}}, \code{\link{chipseq}}, \code{\link{ghist}},
##' \code{\link{gdist}}, \code{\link{css}}, \code{\link{elbow}}, \code{\link{heatmap.3}}
##' @examples
##' \dontrun{
##' require(GMD)        # load GMD
##' help(GMD)           # a help document of GMD 
##' data(package="GMD") # a list of datasets available in GMD
##' ls("package:GMD")   # a list of functions available in GMD
##' help(package="GMD") # help documentation on GMD
##' citation("GMD")     # citation for publications
##' demo("GMD-demo")    # run the demo
##'
##' ## view GMD Description
##' packageDescription("GMD")
##' 
##' ## view GMD vignette
##' vignette("GMD-vignette",package="GMD")
##' }
NULL
