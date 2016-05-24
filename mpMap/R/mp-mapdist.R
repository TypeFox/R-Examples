#' Conversion between recombination fractions (R) and map distance (X) 
#'
#' Haldane and Kosambi map functions and inverses for converting between 
#' recombination fractions and map distance (in cM)
#' @export
#' @rdname mp-mapdist
#' @aliases haldaneX2R haldaneR2X kosambiR2X kosambiX2R mp.mapdist
#' @usage haldaneR2X(r)
#' haldaneX2R(x)
#' kosambiR2X(r)
#' kosambiX2R(x)
#' @param x Map distance (measured in centiMorgans)
#' @param r Recombination fraction
#' @references Jurg Ott, Analysis of Human Genetic Linkage, Johns Hopkins University Press, Baltimore 1999

haldaneR2X <-
function(r) return(-50*log(1-2*r))

haldaneX2R <-
function(x) return(.5*(1-exp(-2*x/100)))

kosambiR2X <-
function(r) return(50*atanh(2*r))

kosambiX2R <-
function(x) return(.5*tanh(2*x/100))

