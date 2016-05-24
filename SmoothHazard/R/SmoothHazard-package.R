#' 
#' Paquid data set
#' 
#' Paquid data set composed of 1000 subjects selected randomly from the Paquid
#' data set of 3675 subjects.
#' 
#' 
#' @name Paq1000
#' @docType data
#' @format A data frame with 1000 rows and the following 8 columns.  \describe{
#' \item{dementia}{dementia status, 0=non-demented, 1=demented}
#' \item{death}{death status, 0=alive, 1=dead} \item{e}{age at
#' entry in the study} \item{l}{for demented subjects: age at the visit
#' before the diagnostic visit; for non-demented subjects: age at the last
#' visit (censoring age)} \item{r}{for demented subjects: age at the
#' diagnostic visit; for non-demented subjects: age at the last visit
#' (censoring age)} \item{t}{for dead subjects: age at death; for alive
#' subject: age at the latest news} \item{certif}{primary school
#' certificate:\code{0=with certificate}, \code{1=without certificate}}
#' \item{gender}{gender: \code{0=female}, \code{1=male}} }
#' @keywords datasets
#' @examples
#' 
#' data(Paq1000)
#' 
NULL



#' Data set for survival models: right-censored and interval-censored data.
#' 
#' A simulated data frame for survival models composed of right-censored and
#' interval-censored data.
#' 
#' 
#' @name testdata
#' @docType data
#' @format A data frame with 936 observations on the following 4 variables.
#' \describe{ \item{l}{for diseased subjects: left endpoint of
#' censoring interval; for non-diseased subjects: right censoring time}
#' \item{r}{for diseased subjects: right endpoint of censoring
#' interval; for non-diseased subjects: right censoring time for the disease
#' event} \item{id}{disease status} \item{cov}{covariate} }
#' @keywords datasets
#' @examples
#' 
#' data(testdata)
#' head(testdata)
#' 
NULL


#' @importFrom lava sim regression
