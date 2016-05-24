#' CDS function definition
#'
#' Creates model formula list for conventional distance sampling using values
#' supplied in call to \code{\link{ddf}}
#'
#' @param key string identifying key function (currently either "hn"
#'   (half-normal),"hr" (hazard-rate), "unif" (uniform) or "gamma" (gamma
#'   distribution)
#' @param adj.series string identifying adjustment functions cos (Cosine), herm
#'   (Hermite polynomials), poly (simple polynomials) or NULL
#' @param adj.order vector of order of adjustment terms to include
#' @param adj.scale whether to scale the adjustment terms by "width" or "scale"
#' @param adj.exp if TRUE uses exp(adj) for adjustment to keep f(x)>0
#' @param formula formula for scale function (included for completeness only
#'   only formula=~1 for cds)
#' @param shape.formula formula for shape function
#' @return A formula list used to define the detection function model
#'   \item{fct}{string "cds"} \item{key}{key function string}
#'   \item{adj.series}{adjustment function string} \item{adj.order}{adjustment
#'   function orders} \item{adj.scale}{adjustment function scale type}
#'   \item{formula}{formula for scale function} \item{shape.formula}{formula
#'   for shape function}
#' @author Jeff Laake; Dave Miller
#' @keywords utility
cds <- function(key=NULL,adj.series=NULL,adj.order=NULL,adj.scale="width",
                adj.exp=FALSE,formula=~1,shape.formula=~1){
# Since we only have a special case of mcds here, lets just let it
# do the work.
  return(mcds(formula,shape.formula=shape.formula,key=key,adj.series=adj.series,adj.order=adj.order,adj.scale=adj.scale,adj.exp=adj.exp))
}
