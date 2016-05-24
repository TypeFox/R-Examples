#' Autism data
#' 
#' Data from a prospective longitudinal study following 214 children between 
#' the ages of 2 and 13 who were diagnosed with either autism spectrum disorder
#' or non-spectrum developmental delays at age 2.
#' 
#' @usage data(autism)
#' @format A data frame with 604 observation on the following 7 variables:
#' \describe{
#' \item{childid}{Child ID.}
#' \item{sicdegp}{Sequenced Inventory of Communication Development group (an
#'   assessment of expressive language development) - a factor. Levels are 
#'   \code{low}, \code{med}, and \code{high}.}
#' \item{age2}{Age (in years) centered around age 2 (age at diagnosis).}
#' \item{vsae}{Vineland Socialization Age Equivalent}
#' \item{gender}{Child's gender - a factor. Levels are \code{male} and \code{female}.}
#' \item{race}{Child's race - a factor. Levels are \code{white} and \code{nonwhite}.}
#' \item{bestest2}{Diagnosis at age 2 - a factor. Levels are \code{autism} and 
#'   \code{pdd} (pervasive developmental disorder).}
#' }
#' @name autism
#' @docType data
#' @keywords datasets
#' @source \url{http://www-personal.umich.edu/~kwelch/}
#' @references 
#'   Anderson, D. K., Lord, C., Risi, S., DiLavore, P. S., Shulman, C., 
#'   Thurm, A., et al. (2007). Patterns of growth in verbal abilities among 
#'   children with autism spectrum disorder. \emph{Journal of Consulting and 
#'   Clinical Psychology}, \bold{75}(4), 594--604. 
#'   
#'   Anderson, D. K., Oti, R. S., Lord, C., & Welch, K. (2009). 
#'   Patterns of Growth in Adaptive Social Abilities Among Children with Autism 
#'   Spectrum Disorders. \emph{Journal of Abnormal Child Psychology}, 
#'   \bold{37}(7), 1019--1034.
NULL