#' Community Dynamics Metrics
#' @description Temporal diversity indices and community stability metrics for ecologists.
#' @details The functions in \code{codyn} implement metrics that are explicitly temporal, and include the option to  calculate them over multiple replicates.
#' Functions fall into two categories: temporal diversity indices and community stability metrics.
#' The diversity indices in \code{codyn} are temporal analogs to traditional diversity indices such as richness and rank-abundance curves.
#' Specifically, \code{codyn} includes functions to calculate species turnover, mean rank shifts and lags in community similarity between time points.
#' The community stability metrics in \code{codyn} calculate overall stability and patterns of species covariance and synchrony over time.
#' Finally, \code{codyn} contains vignettes that describe methods and reproduce figures from published papers to help users contextualize and apply functions to their own data.
#' Work on this package was supported by NSF-ABI grant #1262458
#'
#' @author
#' \itemize{
#'  \item{Lauren Hallett \email{lauren.m.hallett@@gmail.com}}
#'  \item{Sydney K. Jones \email{syd@@sevilleta.unm.edu}}
#'  \item{A. Andrew A. MacDonald \email{aammacdonald@@gmail.com}}
#'  \item{Matthew B. Jones \email{jones@@nceas.ucsb.edu}}
#'  \item{Dan F. B. Flynn \email{flynn@@fas.harvard.edu}}
#'  \item{Peter Slaughter \email{slaughter@@nceas.ucsb.edu}}
#'  \item{Corinna Gries \email{cgries@@wisc.edu}}
#'  \item{Scott L. Collins \email{scollins@@sevilleta.unm.edu}}
#' }
#'
#' @docType package
#' @name codyn
#' @aliases codyn
#' @section Functions:
#' \itemize{
#'  \item{\code{\link[=turnover]{turnover}}}{: Calculates species turnover between time periods}
#'  \item{\code{\link[=rank_shift]{rank_shift}}}{: Calculates the mean relative change in species rank abundances}
#'  \item{\code{\link[=rate_change]{rate_change}}}{: Calculates the rate change in a community over time}
#'  \item{\code{\link[=rate_change_interval]{rate_change_interval}}}{: Produces a data frame containing differences in species composition between samples at increasing time intervals}
#'  \item{\code{\link[=community_stability]{community_stability}}}{: Calculates community stability over time}
#'  \item{\code{\link[=variance_ratio]{variance_ratio}}}{: Computes the ratio of the variance of aggregate species abundances in a community }
#'  \item{\code{\link[=synchrony]{synchrony}}}{: Calculates the degree synchrony in species abundances}
#'  \item{\code{\link[=cyclic_shift]{temporal_torus_translation}}}{: Calculates a test statistic on a null ecological community created via cyclic shifts. \code{confint} provides mean and confidence intervals of this null distribution}
#' }
NULL

#' Konza data from Collins et al. 2008
#'
#' A dataset of tallgrass prairie plant composition at one annually burned and one unburned
#' site over time at the Konza Prairie LTER, Manhattan Kansas (Collins et al. 2008).
#'
#' A data frame containing a column for replicate, year, species and abundance:
#' \itemize{
#'   \item replicate: A factor column of spatial replicates with two levels ("annually burned" and "unburned")
#'   \item year: An integer column of sampling time points
#'   \item species: A factor column of species sampled
#'   \item abundance: A numeric column of abundance values
#' }
#'
#' @source
#' Collins, Scott L., Katharine N. Suding, Elsa E. Cleland, Michael Batty, Steven C. Pennings, Katherine L. Gross, James B. Grace, Laura Gough, Joe E. Fargione, and Christopher M. Clark. (2008) "Rank clocks and plant community dynamics." Ecology 89, no. 12: 3534-41.
#' @docType data
#' @keywords datasets
#' @name collins08
#' @usage data(collins08)
#' @format A data frame with 2058 rows and 4 variables
NULL

#' Data from Konza Prairie, watershed 001d
#'
#' Plant composition within multiple replicates at an annually burned tallgrass
#' prairie site in the Konza Prairie LTER, Manhattan KS (Watershed 001d).
#'
#' A data frame containing a column for species, year, subplot and abundance:
#' \itemize{
#'   \item species: A factor column of species sampled
#'   \item year: An integer column of sampling time points
#'   \item subplot: A factor column of spatial replicates with 20 levels
#'   \item abundance: A numeric column of abundance values
#' }
#'
#' @source
#' Konza Prairie LTER Dataset ID: PVC02, watershed 1D
#'
#' Collins, S. L. (2000) Disturbance frequency and community stability in native tallgrass prairie. American Naturalist 155:311-325.
#' @docType data
#' @keywords datasets
#' @name knz_001d
#' @usage data(knz_001d)
#' @format A data frame with 8768 rows and 4 variables
NULL