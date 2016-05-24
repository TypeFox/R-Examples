#' @name dcars
#' @title Dutch Cars Data
#' @description This data set relates to 187 Dutch households rating 10 automobile manufacturers according 
#' to 8 variables (original Dutch terms in parentheses): price (prijsniveau), design (vormgeving), safety (veiligheid), operating cost (gebruikskosten), )
#' sportiness (sportiviteit), size (modelgrootte), reliability (betrouwbaarheid) and feautures (uitrusting). 
#' A rating scale from 1 to 10 was used. 
#' @details The original sample consisted of 188 households. However, one of these households (code 87845) was discarded 
#' because it appears that they used a rating scale from 0 to 10 instead of from 1 to 10. Note that all rating scales has
#' been reversed so that higher scores are better for most items. The exceptions are OperatingCost and Size, where larger
#' values mean higher costs and smaller cars respectively.
#' @docType data
#' @usage dcars
#' @format A three-way array with cars in the first dimension, variables in the second and 
#' consumers in the third dimension.
#' \describe{
#' The items and labels for the endpoints of the scales are (original Dutch labels in parentheses):
#'  \item{Affordability}{A rating from 1 = Expensive (duur) to 10 = Cheap (goedkoop)}
#'  \item{Attractiveness}{A rating from 1 = Ugly (lelijk) to 10 = Beautiful (mooi)}
#'  \item{Safety}{A rating from 1 = Bad (slecht) to 10 = Good (goed)}
#'  \item{OperatingCost}{A rating from 1 = Low (laag) to 10 = High (hoog)}
#'  \item{Sportiness}{A rating from 1 = Slow (langzaam) to 10 = Fast (snel)}
#'  \item{Size}{A rating from 1 = Large (groot) to 10 = Small (klein)}
#'  \item{Reliability}{A rating from 1 = Bad (slecht) to 10 = Good (goed)}
#'  \item{Features}{A rating from 1 = Simple (eenvoudig) to 10 = Luxurious (luxe)}
#' }
#' @examples
#' data("dcars")
#' set.seed(5448)
#' m <- lsbclust(data = dcars, delta = c(1, 1, 1, 1), nclust = c(5, 3, 6, 8), nstart = 5, 
#'               nstart.kmeans = 10, parallelize = FALSE, fixed = "columns")
#' @source Tammo Bijmolt, Michel van de Velden
NULL