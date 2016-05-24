#' @title Catch-at-age for Gulf Menhaden, 1964-2004.
#' 
#' @description Estimated catch-at-age for Gulf Menhaden (\emph{Brevoortia patronus}), 1964-2004 from thereduction fishery in the U.S. Gulf of Mexico.
#' 
#' @name Menhaden1
#' 
#' @docType data
#' 
#' @format A data frame with 41 observations on the following 7 variables.
#'  \describe{
#'    \item{year}{Year of capture.}
#'    \item{age0}{Estimated catch (millions) of age-0 fish.}
#'    \item{age1}{Estimated catch (millions) of age-1 fish.}
#'    \item{age2}{Estimated catch (millions) of age-2 fish.}
#'    \item{age3}{Estimated catch (millions) of age-3 fish.}
#'    \item{age4}{Estimated catch (millions) of age-4 fish.}
#'    \item{age5}{Estimated catch (millions) of age-5 fish.}
#'    \item{age6}{Estimated catch (millions) of age-6 fish.}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Mortality
#'    \item Catch curve
#'  }
#'  
#' @concept Mortality 'Catch Curve'
#'
#'  @source From Table 2 in Vaughan, D.S., K.W. Shertzer, and J.W. Smith.  2007.  Gulf menhaden (\emph{Brevoortia patronus}) in the U.S. Gulf of Mexico: Fishery characteristics and biological reference points for management.  Fisheries Research 83:263-275.  [Was (is?) from http://menhaden.gsmfc.org/FishRes_Vaughan_etal_2007-GM.pdf.]
#'  
#' @keywords datasets
#' 
#' @examples
#' data(Menhaden1)
#' str(Menhaden1)
#' head(Menhaden1)
#' ages <- 0:6
#' # Extract one year, delete year column (the -1), and transpose to be a vector
#' ct <- t(Menhaden1[Menhaden1$year==1974,-1])
#' plot(ct~ages,pch=16,type="b",xlab="Age",ylab="Est. Catch (Millions)",main="year==1974")
#' 
NULL