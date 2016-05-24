#' Bordeaux Wines Dataset
#' 
#' Quality measures of wines from Bordeaux, France
#' 
#' 
#' @name bordeaux
#' @docType data
#' @format A data frame with 34 observations on the following 6 variables.
#' \tabular{ll}{ \code{year} \tab year of harvest\cr \code{temperature} \tab
#' sum of daily average temperatures (in celsius degrees)\cr \code{sun} \tab
#' duration of insolation (in hours)\cr \code{heat} \tab number of super-hot
#' days\cr \code{rain} \tab rain level (in millimeters)\cr \code{quality} \tab
#' wine quality: a factor with levels \code{bad}, \code{good}, and
#' \code{medium}\cr }
#' @references Chapter 10: Analyse Discriminante, page 353. \cr Tenenhaus M.
#' (2007) \emph{Statistique}. Dunod, Paris.
#' @keywords datasets
#' @examples
#' 
#'  \dontrun{
#'   # load data
#'   data(bordeaux)
#'   
#'   # structure of data
#'   str(bordeaux)
#'  }
#' 
NULL





#' Tools of the Trade for Discriminant Analysis
#' 
#' DiscriMiner contains several functions for Discriminant Analysis and
#' Classification purposes covering various methods such as descriptive,
#' geometric, linear, quadratic, PLS, as well as qualitative discriminant
#' analyses.
#' 
#' \tabular{ll}{ Package: \tab DiscriMiner\cr Type: \tab Package\cr Version:
#' \tab 0.1-23\cr Date: \tab 2012-12-20\cr License: \tab GPL-3\cr }
#' 
#' @name DiscriMiner-package
#' @docType package
#' @author Gaston Sanchez
#' 
#' Maintainer: Gaston Sanchez <gaston.stat@@gmail.com>
#' @references \url{http://www.gastonsanchez.com/discriminer}
#' 
#' Lebart L., Piron M., Morineau A. (2006) \emph{Statistique exploratoire
#' multidimensionnelle}. Dunod, Paris.
#' 
#' Nakache J-P., Confais J. (2003) \emph{Statistique explicative appliquee}.
#' Editions Technip, Paris.
#' 
#' Saporta G. (2006) \emph{Probabilites, analyse des donnees et statistique}.
#' Editions Technip, Paris.
#' 
#' Tenenhaus M. (1998) \emph{La Regression PLS}. Editions Technip, Paris.
#' 
#' Tenenhaus M. (2007) \emph{Statistique}. Dunod, Paris.
#' 
#' Tuffery S. (2008) \emph{Data Mining et Statistique Decisionnelle}. Editions
#' Technip, Paris.
#' 
#' Tuffery S. (2011) \emph{Data Mining and Statistics for Decision Making}.
#' Wiley, Chichester.
#' 
#' \emph{Multiple Correspondence Analysis and Related Methods}. (2006) Edited
#' by Michael Greenacre and Jorg Blasius. Chapman and Hall/CRC
#' @keywords package
NULL





#' Infarctus dataset
#' 
#' Infarctus dataset from Saporta (2006)
#' 
#' 
#' @name infarctus
#' @docType data
#' @format A data frame with 101 observations on the following 8 variables.
#' \tabular{ll}{ \code{FRCAR} \tab Frequence Cardiaque (i.e. heart rate)\cr
#' \code{INCAR} \tab Index Cardique (cardiac index)\cr \code{INSYS} \tab Index
#' Systolique (systolic index)\cr \code{PRDIA} \tab Pression Diastolique
#' (diastolic pressure)\cr \code{PAPUL} \tab Pression Arterielle Pulmonaire
#' (pulmonary artery pressure)\cr \code{PVENT} \tab Pression Ventriculaire
#' (ventricular pressure)\cr \code{REPUL} \tab Resistance Pulmonaire (pulmonary
#' resistance)\cr \code{PRONO} \tab Pronostic (prognosis): a factor with levels
#' \code{dead} and \code{survive}\cr }
#' @references Chapter 18: Analyse discriminante et regression logistique, pp
#' 453-454 \cr Saporta G. (2006) \emph{Probabilites, analyse des donnees et
#' statistique}. Editions Technip, Paris.
#' @keywords datasets
#' @examples
#' 
#'  \dontrun{
#'   # load data
#'   data(infarctus)
#'   
#'   # summary of variables
#'   summary(infarctus)
#'  }
#' 
NULL





#' Insurance Dataset
#' 
#' Dataset of car-insurance customers from Belgium in 1992
#' 
#' Dataset for DISQUAL method
#' 
#' @name insurance
#' @docType data
#' @format
#' 
#' A data frame with 1106 observations on the following 10 variables.
#' \tabular{ll}{ \code{Claims} \tab Group variable. A factor with levels
#' \code{bad} and \code{good}\cr \code{Use} \tab Type of Use. A factor with
#' levels \code{private} and \code{professional}\cr \code{Type} \tab Insurance
#' Type. A factor with levels \code{companies}, \code{female}, and
#' \code{male}\cr \code{Language} \tab Language. A factor with levels
#' \code{flemish} and \code{french}\cr \code{BirthCohort} \tab Birth Cohort. A
#' factor with levels \code{BD_1890_1949}, \code{BD_1950_1973}, and
#' \code{BD_unknown}\cr \code{Region} \tab Geographic Region. A factor with
#' levels \code{Brussels} and \code{Other_regions}\cr \code{BonusMalus} \tab
#' Level of bonus-malus. A factor with levels \code{BM_minus} and
#' \code{BM_plus}\cr \code{YearSuscrip} \tab Year of Subscription. A factor
#' with levels \code{YS<86} and \code{YS>=86}\cr \code{Horsepower} \tab
#' Horsepower. A factor with levels \code{HP<=39} and \code{HP>=40}\cr
#' \code{YearConstruc} \tab Year of vehicle construction. A factor with levels
#' \code{YC_33_89} and \code{YC_90_91}\cr }
#' @seealso \code{\link{disqual}}
#' @references Saporta G., Niang N. (2006) Correspondence Analysis and
#' Classification. In \emph{Multiple Correspondence Analysis and Related
#' Methods}, M. Greenacre and J. Blasius, Eds., pp 371-392. Chapman & Hall/CRC,
#' Boca Raton, Florida, USA.
#' @keywords datasets
#' @examples
#' 
#'  \dontrun{
#'   # load data
#'   data(insurance)
#'   
#'   # structure
#'   str(insurance)
#'  }
#' 
NULL



