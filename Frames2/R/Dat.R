#' @name Dat
#' @aliases Dat
#' @docType data
#' @title Joint sample database
#' 
#' @description This dataset contains some variables coming from a real dual frame survey conducted in 2013 in Andalusia (Spain) by a scientific institute specialized in social topics.
#' With this dataset it is intented to show how to properly split a joint dual frame sample into subsamples, so functions of \code{Frame2} can be used.  
#' @usage Dat
#' @details The survey was based on two frames: a landline frame and a cell phone frame. Landline frame was stratified by province and simple random sampling without replacement was considered
#' in cell phone frame. The size of the whole sample was \eqn{n = 2402}. Total of the variable Income in the whole population is \eqn{X_{Income} = 12686232063}.
#' @format
#' \describe{
#'    \item{Drawnby}{Indicates whether individual was selected in the landline sample(1) or in the cell phone sample(2).}
#'    \item{Stratum}{Indicates the stratum each individual belongs to. For individuals selected in cell phone sample, value of this variable is \code{NA}.}
#'    \item{Opinion}{Response of the individual to the question: Do you think that immigrants currently living in Andalusia are quite a lot? 1 represents "yes" and 0 represents "no".}
#'    \item{Landline}{Indicates whether individual has a landline (1) or not (0).}
#'    \item{Cell}{Indicates whether individual has a cell phone(1) or not(0).}
#'    \item{ProbLandline}{First order inclusion probability of reaching the individual by landline.}
#'    \item{ProbCell}{First order inclusion probability of reaching the individual by cell phone.}
#'    \item{Income}{Monthly income (in euros) of the individual.}
#' }
#' @examples
#' data(Dat)
#' attach(Dat)
#' 
#' #We are going to split dataset Dat into two new datasets, each 
#' #one corresponding to a frame: frame containing individuals
#' #using landline and frame containing individuals using cell phone.
#' 
#' FrameLandline <- Dat[Landline == 1,]
#' FrameCell <- Dat[Cell == 1,]
#'
#' #Equally, we can split the original dataset in three new different 
#' #datasets, each one corresponding to one domain: first domain containing
#' #individuals using only landline, second domain containing individuals
#' #using only cell phone and the third domain containing individuals
#' #using both landline and cell phone.
#' 
#' DomainLandline <- Dat[Landline == 1 & Cell == 0,]
#' DomainCell <- Dat[Landline == 0 & Cell == 1,]
#' DomainBoth <- Dat[Landline == 1 & Cell == 1,]
#'
#' #From the domain datasets, we can build frame datasets
#'
#' FrameLandline <- rbind(DomainLandline, DomainBoth)
#' FrameCell <- rbind(DomainCell, DomainBoth)
#' @keywords datasets
NULL