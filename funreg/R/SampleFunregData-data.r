#' @name SampleFunregData
#' @title Sample dataset for funreg
#' @description A data frame generated using the following code:
#'  \code{set.seed(123);
#'        SampleFunregData <- generate.data.for.demonstration(nsub = 200,  
#'                                                            b0.true = -5,
#'                                                            b1.true = 0, 
#'                                                            b2.true = +1,
#'                                                            b3.true = -1,
#'                                                            b4.true = +1,
#'                                                            nobs = 400, 
#'                                                            observe.rate = 0.1);}
#' @docType data
#' @format A data frame for a simulated longitudinal study, 
#' in "tall" rather than "wide" format (multiple rows per
#' individual, one for each measurement time)
#' with 8109 rows and 13 columns.
#'
#' \describe{
#' 
#' \item{id}{Integer uniquely identifying the subject to whom this data row pertains.}
#' 
#' \item{s1,s2,s3,s4}{Four subject-level (time-invariant) covariates.}
#' 
#' \item{y}{A response, coded as 0 or 1, which is to be modeled 
#'  using a functional regression.  It is also subject-level (i.e., 
#'  either time-invariant or measured only once).  However, 
#'  like \code{s1} through \code{s4}, its value is repeated for each 
#'  row of data for a subject.}
#' 
#' \item{time}{The time variable, arbitrarily chosen to range
#'  from a low of 0 to a high of 10, which identifies when this row's 
#'  observations are taken.}
#' 
#' \item{true.x1,true.x2}{The unknown smooth expected values of two
#' time-varying variables which can be treated as functional covariates. They
#' vary by subject and time and are therefore different in each row.}
#' 
#' \item{true.betafn1,true.betafn2}{The unknown true functional
#' regression coefficient function used to generate \code{y} from the two time-varying 
#' predictors.  The latter is always zero because x2 is unrelated to y.}
#' 
#' \item{x1,x2}{The observed values of the two functional regression
#' predictors, as measured for a given time on a given subject.}
#' }
#' @keywords datasets 
#' 
 NULL  # see http://stackoverflow.com/questions/7086805/how-can-i-document-datasets-without-adding-them-to-the-collate-field