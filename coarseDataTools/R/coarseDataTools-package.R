##' @name exp.win.lengths
##' @title Exposure window lengths from an influenza outbreak at a NYC school
##' @description A numeric vector of exposure window lengths taken from a dataset of doubly interval-censored incubation period observations.  All observations came from a NYC public school.  The outbreak has been described in full in Lessler et al. (see citation below). 
##' @docType data
##' @format A numeric vector with 134 positive values.  Each value represents an exposure window length from an observation of the incubation period for that individual.  The exposure window length is the length of time during which exposure could have occured.  For example, if an individual could have been exposed anytime between 6am on Monday to 6am on Wednesday, her exposure window length would be 2 days. 
##' 
##' @usage data(exp.win.lengths)
##' @source Lessler J et al.  New England Journal of Medicine. Outbreak of 2009 Pandemic Influenza A (H1N1) at a New York City School. 2009. 361(27):2628-2636. \url{http://content.nejm.org/cgi/content/full/361/27/2628}
##' 
##' @examples
##' data(exp.win.lengths)
##' summary(exp.win.lengths)
##' hist(exp.win.lengths)
##' @keywords datasets
NULL

##' @name fluA.inc.per
##' @docType data
##' @title Coarse incubation period data for influenza A
##' @description These observations on the incubation period of influenza A come from a variety of sources, and were gathered for a literature review.  They report doubly interval-censored, single interval-censored or exact observations for the incubation period. 
##' @usage data(fluA.inc.per)
##' @format A data frame with 151 observations on the following 7 variables. 
##' \describe{
##' \item{\code{author}}{the name of the primary author for the source of the observation}
##' \item{\code{year}}{the year of the study which is the source of the observation}
##' \item{\code{EL}}{the earliest possible time of infection}
##' \item{\code{ER}}{the latest possible time of infection}
##' \item{\code{SL}}{the earliest possible time of symptom onset}
##' \item{\code{SR}}{the latest possible time of symptom onset}
##' \item{\code{type}}{an indicator of the type of observation: 0 for doubly interval-censored, 1 for single-interval censored, 2 for exact}}
##'@source Lessler J, Reich NG, Brookmeyer R, Perl TM, Nelson KE, Cummings DAT. (2009) A systematic review of the incubation periods of acute respiratory viral infections. Lancet Infectious Diseases. 9(5):291-300.
##' @examples
##' data(fluA.inc.per)
##' head(fluA.inc.per)
##' @keywords datasets
NULL

##' @name simulated.outbreak.deaths
##' @docType data
##' @title Simulated case and death reports from a fictional outbreak
##' @description This dataset provides reported counts of cases and deaths occuring at different time points across a simulated outbreak. Details of the data simulation algorithm are provided in the manuscript "Estimating case fatality ratios from infectious disease surveillance data" (Reich et al., under review, available upon request).
##' @usage data(simulated.outbreak.deaths)
##' @format
##' \describe{
##'             \item{\code{time}}{time, t,  after start of outbreak}
##'             \item{\code{grp}}{an categorical variable indicating membership in one of two groups of a covariate, j}
##'             \item{\code{R}}{number of recovered cases reported at the given t and j}
##'             \item{\code{D}}{number of deaths reported at the given t and j}
##'             \item{\code{N}}{total number of cases and deaths reported at t and j, or D+R}
##'             }
##' @source Reich NG, Lessler J, Cummings DAT, Brookmeyer R. Estimating case fatality ratios from infectious disease surveillance data. Biometrics. 2012. 68(2): 598-606.
##' @examples
##' data(simulated.outbreak.deaths)
##' head(simulated.outbreak.deaths)
##' plot(simulated.outbreak.deaths[simulated.outbreak.deaths[,"grp"]==1,"D"], type="l")
NULL

#' @name nycH1N1
##' @docType data
##' @title Incubation period data from New York City Public Schools, 2009 H1N1 influenza outbreak
##' @description These observations on the incubation period of influenza A come from the investigation of the H1N1 outbreak in NYC schools in the spring of 2009. They report doubly interval-censored observations for the incubation period. 
##' @usage data(nycH1N1)
##' @format A data frame with 134 observations on the following 5 variables. 
##' \describe{
##' \item{\code{EL}}{the earliest possible time of infection}
##' \item{\code{ER}}{the latest possible time of infection}
##' \item{\code{SL}}{the earliest possible time of symptom onset}
##' \item{\code{SR}}{the latest possible time of symptom onset}
##' \item{\code{type}}{an indicator of the type of observation: 0 for doubly interval-censored, 1 for single-interval censored, 2 for exact. All of these observations are doubly interval-censored.}}
##' @source Lessler J, Reich NG, Cummings DAT and The DOHMH Swine Influenza Investigation Team. Outbreak of 2009 Pandemic Influenza A (H1N1) at a New York City School. New England Journal of Medicine. 2009. 361(27):2628-2636.
##' @examples
##' data(nycH1N1)
##' head(nycH1N1)
##' dic.fit(nycH1N1)
##' @keywords datasets
NULL