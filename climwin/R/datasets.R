#' Body mass data since 1979.
#' 
#' Artificially generated data representing
#' average body mass since 1979.
#' @format A data frame with 47 rows and 2 variables
#' \describe{
#'   \item{Date}{Date of mass measurements (dd/mm/yyyy).}
#'   \item{Mass}{Annual average body mass in grams.}
#'   }
#'@name Mass
NULL

#' Daily climate data since 1979.
#' 
#' Daily temperature and rainfall data since 1979.
#' @format A data frame with 17,532 rows and 3 variables.
#' \describe{
#'   \item{Date}{Date when climate data was recorded (dd/mm/yyyy).}
#'   \item{Rain}{Daily rainfall data in mm.}
#'   \item{Temp}{Daily temperature data in degrees centigrade.}
#'   }
#'@name MassClimate
NULL

#' Example output dataframe from function climatewin.
#' 
#' Output file from \code{\link{climatewin}} using temperature and body mass data. 
#' Generated with \code{\link{Mass}} and \code{\link{MassClimate}} dataframes.
#' @format A data frame with 5,151 rows and 17 variables.
#' \describe{
#'   \item{deltaAICc}{Difference between model AICc of fitted climate window and a null model containing no climate.}
#'   \item{WindowOpen}{The start day of each tested climate window. Furthest from the biological record.}
#'   \item{WindowClose}{The end day of each tested climate window. Closest to the biological record.}
#'   \item{ModelBeta}{Beta estimate of the relationship between temperature and mass.}
#'   \item{ModelBetaQ}{Quadratic beta estimate of the relationship between temperature and mass.}
#'   \item{ModelBetaC}{Cubic beta estimate of the relationship between temperature and mass.}
#'   \item{ModelInt}{Model intercept.}     
#'   \item{Function}{The function used to fit climate (e.g. linear ("lin"), quadratic ("quad"))}
#'   \item{furthest}{Furthest day back considered in climatewin.}
#'   \item{closest}{Closest day back considered in climatewin.}
#'   \item{Statistics}{The aggregate statistic used to analyse climate (e.g. mean, max, slope).}
#'   \item{type}{Whether "fixed" or "variable" climate windows were tested.}
#'   \item{K}{Number of folds used for k-fold cross validation.}
#'   \item{ModWeight}{Model weight of each fitted climate window.}
#'   \item{cutoff.day,cutoff.month}{If type is "fixed", the date from which the climate window was tested.}
#'   \item{Randomised}{Whether the data was generated using \code{\link{climatewin}} or \code{\link{randwin}}.}
#'   }
#'@name MassOutput
NULL


#' Example output dataframe from function randwin.
#' 
#' Output file from function \code{\link{randwin}} using temperature and mass data. 
#' Generated with \code{\link{Mass}} and \code{\link{MassClimate}} dataframes.
#' @format A data frame with 25,755 rows and 18 variables.
#' \describe{
#'   \item{deltaAICc}{Difference between model AICc of fitted climate window and a null model containing no climate.}
#'   \item{WindowOpen}{The start day of each tested climate window. Furthest from the biological record.}
#'   \item{WindowClose}{The end day of each tested climate window. Closest to the biological record.}
#'   \item{ModelBeta}{Beta estimate of the relationship between temperature and mass.}
#'   \item{ModelBetaQ}{Quadratic beta estimate of the relationship between temperature and mass.}
#'   \item{ModelBetaC}{Cubic beta estimate of the relationship between temperature and mass.}
#'   \item{ModelInt}{Model intercept.}     
#'   \item{Function}{The function used to fit climate (e.g. linear ("lin"), quadratic ("quad"))}
#'   \item{furthest}{Furthest day back considered in climatewin.}
#'   \item{closest}{Closest day back considered in climatewin.}
#'   \item{Statistics}{The aggregate statistic used to analyse climate (e.g. mean, max, slope).}
#'   \item{type}{Whether "fixed" or "variable" climate windows were tested.}
#'   \item{K}{Number of folds used for k-fold cross validation.}
#'   \item{ModWeight}{Model weight of each fitted climate window.}
#'   \item{cutoff.day,cutoff.month}{If type is "fixed", the date from which the climate window was tested.}
#'   \item{Randomised}{Whether the data was generated using \code{\link{climatewin}} or \code{\link{randwin}}.}
#'   \item{Repeat}{The number of randomisations carried out.}
#'   }
#'@name MassRand
NULL
   
#' Reproductive success of birds since 2009.
#' 
#' Artificially generated data representing
#' reproductive success of birds since 2009.
#' @format A data frame with 1,619 rows and 4 variables.
#' \describe{
#'   \item{Offspring}{Total number of offspring produced.}
#'   \item{Date}{Date of hatching (dd/mm/yyyy).}
#'   \item{Order}{Order of nest within each season.}
#'   \item{BirdID}{Individual ID of female.}
#'    }
#'@name Offspring
NULL

#' Daily climate data since 2009.
#' 
#' Daily temperature and rainfall data since 2009.
#' Coincides with biological data from \code{\link{Offspring}}.
#' @format A data frame with 2,588 rows and 3 variables.
#' \describe{
#'   \item{Date}{Date when climate was recorded (dd/mm/yyyy).}
#'   \item{Rain}{Daily rainfall data in mm.}
#'   \item{Temperature}{Daily temperature data in degrees centigrade.}
#'   }
#'@name OffspringClimate
NULL