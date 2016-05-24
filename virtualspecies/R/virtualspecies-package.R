#' Generation of virtual species
#'
#' This package allows generating virtual species distributions, for example
#' for testing species distribution modelling protocols. For a complete tutorial, 
#' see http://borisleroy.com/virtualspecies
#' 
#' 
#' @details
#' The process of generating a virtual species distribution is divided into 
#' four major steps.
#' \enumerate{
#' \item{Generate a virtual species distributions from environmental variables.
#' This can be done by 
#' \itemize{
#' \item{defining partial response functions to each environmental 
#' variable, and then combining them to compute the overall environmental suitability,
#' with \code{\link{generateSpFromFun}}}
#' \item{computing a PCA among
#' environmental variables, and simulating the response of the species along 
#' the two first axes of the PCA with \code{\link{generateSpFromPCA}}}}
#' This step can be randomised with \code{\link{generateRandomSp}}
#' }
#' \item{Convert the virtual species distribution into presence-absence, with
#' \code{\link{convertToPA}}}
#' \item{Facultatively, introduce a distribution bias with 
#' \code{\link{limitDistribution}}}
#' \item{Sample occurrence points (presence only or presence-absence) inside the
#' virtual species distribution, either randomly or with biases, with 
#' \code{\link{sampleOccurrences}}}
#' }
#' 
#' There are other useful functions in the package:
#' \itemize{
#' \item{\code{\link{formatFunctions}}: this is a helper function to format and 
#' illustrate the response functions as a correct input for \code{\link{generateSpFromFun}}}
#' \item{\code{\link{plotResponse}}: to visualise the species-environment relationship
#' of the virtual species}
#' \item{\code{\link{removeCollinearity}}: this function can be used to remove
#' collinearity among variables of a stack by selecting a subset of 
#' non-intercorrelated variables}
#' \item{\code{\link{synchroniseNA}}: this function can be used to synchronise
#' NA values among layers of a stack}
#' }
#' 
#' 
#' 
#' 
#' This packages makes use of different other packages:
#' \itemize{
#' \item{This package makes extensive use of the package \code{\link{raster}} to obtain spatialised
#' environmental variables, and produce spatialised virtual species distributions.}
#' \item{\code{\link{dismo}} is used to sample occurrence data from the generated virtual species.}
#' \item{\code{\link{ade4}} is used to generate species with a PCA approach.}
#' \item{\code{\link{rworldmap}} is used to obtain free world shapefiles, in order to create
#' dispersal limitations and sampling biases.}
#' }
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#' 
#' with help from C. N. Meynard, C. Bellard & F. Courchamp
#' 
#' Maintainer: Boris Leroy \email{leroy.boris@@gmail.com}
#' @name virtualspecies-package
#' @docType package
NULL