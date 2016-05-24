#' Database of embryonic development and thermosensitive period of development for sex determination
#' @title Database of of embryonic development and thermosensitive period of development for sex determination
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name stages
#' @description Database of embryonic development and thermosensitive period of development for sex 
#' determination.
#' @references Pieau, C., Dorizzi, M., 1981. Determination of temperature sensitive 
#' stages for sexual differentiation of the gonads in embryos of the turtle, 
#' Emys orbicularis. Journal of Morphology 170, 373-382.\cr
#' Yntema, C.L., Mrosovsky, N., 1982. Critical periods and pivotal temperatures for 
#' sexual differentiation in loggerhead sea turtles. Canadian Journal of 
#' Zoology-Revue Canadienne de Zoologie 60, 1012-1016.\cr
#' Kaska, Y., Downie, R., 1999. Embryological development of sea turtles (Chelonia mydas, 
#' Caretta caretta) in the Mediterranean. Zoology in the Middle East 19, 55-69.\cr
#' Greenbaum, E., 2002. A standardized series of embryonic stages for the emydid 
#' turtle Trachemys scripta. Canadian Journal of Zoology-Revue Canadienne de 
#' Zoologie 80, 1350-1370.
#' @keywords datasets
#' @family Functions for temperature-dependent sex determination
#' @usage stages
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(stages)
#' names(stages)
#' levels(as.factor(stages$Species))
#' }
#' @format A list with dataframes including attributes
NULL
