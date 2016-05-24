#' Replication data for Clark and Golder (2006)
#'
#' District magnitude and ethnic heterogeneity data from a pooled sample of 
#' established democracies in the 1990s. Data originally
#' from Clark and Golder (2006).
#'
#' @format A data frame with 754 rows and 33 variables:
#' \describe{
#'   \item{country}{country name}
#'   \item{countrynumber}{country number}
#'   \item{year}{year of observation}
#'   \item{enep1}{electoral parties}
#'   \item{eneg}{ethnic heterogeneity}
#'   \item{logmag}{district magnitude}
#'   \item{legelec}{legislative election}
#'   \item{preselec}{presidential election}
#'   \item{regime}{regime as of 31 Dec of given year (0=democracy, 1=dictatorship)}
#'   \item{regime_leg}{regime type at time of leg. election (0=democracy, 1=dictatorship)}
#'   \item{eighties}{election in 1980s closest to 1985}
#'   \item{nineties}{election in 1990s closest to 1995}
#'   \item{old}{elections in countries that did not transition to democracy in 1990s}
#'   \item{avemag}{average district magnitude}
#'   \item{districts}{number of electoral districts}
#'   \item{enep}{effective number of ethnic groups fearon}
#'   \item{enep_others}{n/a}
#'   \item{enpp}{parliamentary parties - uncorrected}
#'   \item{enpp_others}{n/a}
#'   \item{enpp1}{parliamentary parties - corrected}
#'   \item{enpres}{effective number of presidential candidates}
#'   \item{medmag}{median district magnitude}
#'   \item{newdem}{first election of new democracy}
#'   \item{proximity1}{proximity - continuous}
#'   \item{proximity2}{proximity - dichotomous}
#'   \item{seats}{assembly size}
#'   \item{upperseats}{number of upper tier seats}
#'   \item{uppertier}{percentage of uppertier seats}
#'   \item{uppertier_eneg}{uppertier*eneg}
#'   \item{logmag_eneg}{logmag*eneg}
#'   \item{proximity1_enpres}{proximity1*enpres}
#'   \item{twoelections}{n/a}
#'   \item{twoelections1}{n/a}
#'   ...
#' }
#' @source Clark, William R., and Matt Golder. 2006. "Rehabilitating Duverger's Theory." \emph{Comparative Political Studies} 39(6): 679-708.
#' @name legfig
NULL