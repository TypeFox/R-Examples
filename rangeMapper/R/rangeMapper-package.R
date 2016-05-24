

#' rangeMapper: A platform for the study of macroecology of life history
#' traits.
#'
#' \pkg{rangeMapper} is a front end platform for the study of macroecology of
#' life history traits at both inter-specific and assemblage levels.
#'
#' @name       rangeMapper
#' @docType    package
#' @references Valcu, M., Dale, J. and Kempenaers, B. (2012) rangeMapper: A
#'      platform for the study of macroecology of life history traits. 21(9). (DOI:
#'      10.1111/j.1466-8238.2011.00739.x)

NULL

#' Life history data of the New World Wrens
#'
#' Life history data (body size, body mass and clutch size) of 84 wren
#' \strong{(Troglodytidae)} species
#'
#' Taxonomic nomenclature follows (Kroodsma & Brewer, 2005) with the exception
#' of \emph{Donacoblis atricapilla} which has been excluded due to its
#' uncertain taxonomic position.
#'
#' @name     wrens
#' @docType  data
#' @format   A data frame with 84 observations on the following 7 variables.
#' \describe{
#'    \item{ID_HBW}{ Handbook of the birds of the world ID}
#'    \item{body_mass}{body mass (grams)}
#'    \item{body_size}{body size (cm)}
#'    \item{clutch_size}{mean or modal clutch size}
#'    \item{com_name}{English name; a factor with 84 levels}
#'    \item{genus}{Genus name}
#'    \item{sci_name}{scientific name; a factor with 84 levels}
#'    \item{source}{bibliographic source of each trait (see references)} }
#' @seealso \code{\link{rangeMap.save}}.
#' @references Auer, S.K., Logue, D.M., Bassar, R.D. & Gammon, D.E. (2007)
#' Nesting biology of the Black-bellied Wren (Thryothorus fasciatoventris) in
#' central Panama. Wilson Journal of Ornithology, 119, 71-76.\cr Dunning, J.B.
#' (2008) CRC handbook of avian body masses, 2nd edn. CRC Press, Boca Raton.\cr
#' Freeman, B.G. & Greeney, H.F. (2008) First description of the nest, eggs and
#' cooperative breeding behavior in sharpe's wren (Cinnycerthia olivascens).
#' Ornitologia Colombiana, 7, 88-92.\cr Hron, K., Templ, M. & Filzmoser, P.
#' (2010) Imputation of missing values for compositional data using classical
#' and robust methods. Computational Statistics & Data Analysis, 54, 3095-3107
#' (function impKNNa using default arguments).\cr Kroodsma, D.E. & Brewer, D.
#' (2005) Family Troglodytidae (Wrens). Lynx Edicions, Barcelona, Spain.\cr
#' Londono, G.A. (2009) Eggs, Nests, and Incubation Behavior of the Moustached
#' Wren (Thryothorus genibarbis) in Manu National Park, Peru. Wilson Journal of
#' Ornithology, 121, 623-627.\cr Ridgely, R.S., T. F. Allnutt, T. Brooks, D. K.
#' McNicol, D. W. Mehlman, B. E. & Young, a.J.R.Z. (2007) Digital Distribution
#' Maps of the Birds of the Western Hemisphere, version 3.0. NatureServe,
#' Arlington, Virginia, USA.\cr Vargas-Soriano, J., Ortiz, J.S. & Segura, G.E.
#' (2010) Breeding Phenology and Nesting Success of the Yucatan Wren in the
#' Yucatan Peninsula, Mexico. Wilson Journal of Ornithology, 122, 439-446.\cr
#' @keywords datasets
#' @examples
#'
#' data(wrens)
#' plot(body_size ~ body_mass, wrens)
#' plot(clutch_size ~ log(body_mass), wrens)
#'
NULL



