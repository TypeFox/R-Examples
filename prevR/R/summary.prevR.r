#' Detailed summary of the variables of a prevR object
#' 
#' Method \code{summary} for objects of class \code{\link[=prevR-class]{prevR}}: 
#' shows a summary of the variables of the object.
#' 
#' @param object object of class \code{\link[=prevR-class]{prevR}}.
#' @param probs vector of probabilities with values in [0,1] for computing 
#'   quantiles of the rings radii (see examples).
#'   
#' @seealso \code{\link{print,prevR-method}}.
#' @examples 
#' summary(fdhs)
#' \dontrun{
#'  dhs <- rings(fdhs, N=c(100,300,500))
#'  summary(dhs)
#'  summary(dhs, c(0, 0.25, 0.5, 0.75, 1))
#' }
#'
#' @aliases summary-methods summary,prevR-method summary prevRsummary

setMethod("summary","prevR",
  function(object, probs = c(0,.10,.25,.50,.75,.80,.90,.95,.99,1)){
    prevRsummary(object,probs)
  }
)

prevRsummary <- function(object, probs = c(0,.10,.25,.50,.75,.80,.90,.95,.99,1)){
  ##################################################################################################
  # Cette fonction renvoie
  # Un resume de l'element clusters (summary.data.frame)
  # Pour chaque ring les quantiles des rayons (pour une probabilite definie par l'argument probs)
  # Un resume statistique de chaque ring (summary.data.frame) 
  ##################################################################################################
  message("Object of class 'prevR'\n",domain="R-prevR")
  clusters = slot(object,"clusters")
  message("SLOT CLUSTERS\n",domain="R-prevR")
  print(summary(clusters[,-1])) # Ne pas afficher la colonne id
  out = NULL
  if(is.prevR(object,"rings")){
    rings = slot(object,"rings")
    for(ring in rings){
      message(gettextf("\nSLOT RINGS FOR N=%s AND R=%s\n",ring$N,ring$R,domain="R-prevR"))
      print(summary(ring$estimates[,-1]))  # Ne pas afficher la colonne id
      
    }
    if (!is.null(probs)){
      projCRS = slot(object,"proj")
      proj    = slot(projCRS,"projargs")
      if(regexpr("longlat",proj)==-1 && regexpr("latlong",proj)==-1){
        r.unit = gettext("in the unit of the projection",domain="R-prevR")
      } else {
        r.unit = gettext("in kilometers",domain="R-prevR")
      }
      message(gettextf("\nQUANTILES OF r.radius (%s):\n",r.unit,domain="R-prevR"))
      sum.rings = function(ring,probs){
        c(ring$N, ring$R, quantile(ring$estimates$r.radius,probs))
      }
      out = t(sapply(rings,sum.rings,probs=probs))
      print(round(out,2)[,c(-1,-2)])
    }
  }
}