#' @name lov
#' @title List-of-values Data Set
#' @description This is the list-of-values data set used in Van Rosmalen, Van Herk & Groenen (2010). 
#' Column names and factor labels differ slightly from that paper. Missing values are encoded as 
#' \code{NA} as usual. The first nine columns are items answered on a nine-point rating scale, with
#' rating 1 representing 'very important' and category 9 'not important at all'. The respondents
#' were asked how important each of these items are as a guiding principle in their lives.
#' @docType data
#' @usage data("lov")
#' @format A data frame with 4514 observations on the following 12 variables.
#' \describe{ 
#'  \item{Belonging}{a numeric vector; 'a sense of belonging'}
#'  \item{Excitement}{a numeric vector} 
#'  \item{Relationships}{a numeric vector; 'warm relationships with others'} 
#'  \item{Self-fulfilment}{a numeric vector}
#' \item{Respected}{a numeric vector; 'being well-respected'} 
#' \item{Enjoyment}{a numeric vector; 'fun and enjoyment'} 
#' \item{Security}{a numeric vector}
#' \item{Self-respect}{a numeric vector}
#' \item{Accomplishment}{a numeric vector; 'a sense of accomplishment'} 
#' \item{Country}{a factor with levels \code{Britain}, \code{France}, \code{Germany}, 
#'    \code{Italy} and \code{Spain}} 
#' \item{Education}{a factor with levels \code{Low} and \code{High}} 
#' \item{Age}{a factor with levels \code{-25}, \code{25-39}, \code{40-54} and \code{55+}}}
#' @source Joost van Rosmalen
#' @references Van Rosmalen, J., Van Herk, H., & Groenen, P. J. (2010). Identifying 
#' response styles: A latent-class bilinear multinomial logit model. 
#' \emph{Journal of Marketing Research}, 47(1), 157-172.
#' @examples
#' data("lov")
#' 
#' ## Construct array
#' lovarr <- indarr(lov[, 1:9], maxcat = 9)
#' 
#' ## Run analysis
#' set.seed(13841)
#' fit <- lsbclust(data = lovarr, margin = 3, delta = c(0, 1, 0, 0), nclust = c(NA, 11, NA, 5), 
#'                  fixed = "rows", nstart = 1, iter.max = 50, nstart.kmeans = 10)
 NULL