#' Test if an object is of class prevR.
#' This function test if the class of an object is \code{\link[=prevR-class]{prevR}}. 
#' It could be used to test the slot \code{rings} or the slot \code{boundary}.
#' 
#' @param object object to test.
#' @param slot "clusters", "rings","boundary" or "proj".
#' @details 
#' Slots \code{rings} and \code{boundary} are always present in an object of class 
#' \code{\link[=prevR-class]{prevR}}, but \code{rings} could be \code{NULL} and 
#' \code{boundary} a \code{\link[sp:SpatialPolygons-class]{SpatialPolygons}} with an 
#' attribute named \code{valid} with the value \code{FALSE} (when boundaries of the studied 
#' area have not been specified explicitly).
#' \itemize{
#'  \item If \code{rings} is \code{NULL}, \code{is.prevR(object,"rings")} will return \code{FALSE}.
#'  \item If \code{boundary} has an attribute \code{valid} equal to \code{FALSE}, 
#'   \code{is.prevR(object,"boundary")} will return \code{FALSE}.
#' } 
#'
#' @return  \code{TRUE} or \code{FALSE}.
#' @seealso \code{\link{prevR-class}}.
#' @examples 
#' col <- c(id = "cluster", 
#'   x = "x",
#'   y = "y",
#'   n = "n",
#'   pos = "pos",
#'   c.type = "residence",
#'   wn = "weighted.n",
#'   wpos = "weighted.pos"
#' )
#' dhs <- as.prevR(fdhs.clusters,col, fdhs.boundary)
#' 
#' is.prevR(dhs)
#' is.prevR(dhs,"rings")
#' is.prevR(dhs,"boundary")
#' 
#' dhs <- rings(dhs,N=300)
#' is.prevR(dhs,"rings")
#'
#' @keywords class
#' @export

is.prevR = function(object, slot = NULL){
  ###############################################################################################
  # Cette fonction permet de tester si l'obejt est de classe prevR
  # Elle permet aussi de tester la presence des elements rings et boundary
  #    Attention les elemeents rings et boundary sont toujours presents mais rings peut etre une liste NULL
  #    et boundary un spatialPolygons avec un attribut valid pose a F (cas ou n'a pas defini de frontieres)
  #    Si rings est une liste NULL is.prevr(object,"rings") renvoie F 
  #    Si boundary a un attribut valid = F is.prevr(object,"boundary") renvoie F 
  # Les arguments 
  #   object : un objet quelconque
  #   slot : un vecteur de chaines de characters contenant clusters, rings, boundary, proj
  #
  # Cette fonction renvoie un vecteur de T ou de F
  #
  # exemple d'utilisation
  # is.prevR(data.prevR)
  # is.prevR(data.prevR,"rings")
  # is.prevR(data.prevR,c("rings","boundary"))
  # 
  ###############################################################################################  
  if(class(object)!="prevR") return(F)
  if(is.null(slot) && class(object)=="prevR") return(T)
  ind = match(slot,slotNames(object),nomatch=0)
  if(any(ind==0)){
    for (i in length(slot[ind==0])) {
      message(gettextf("The slot '%s' doesn't exist for class 'prevR'.",slot[ind==0][i],domain="R-prevR"))
    }
    return(F)
  }
  if(!is.null(slot)){
    response = NULL
    for(one.slot in slot){
      if(one.slot == "boundary"){
        if(attr(slot(object,"boundary"),"valid")){
          response = c(response, T)
        } else {
          response = c(response, F)
        }
        } else {
        if(length(slot(object,one.slot))>0){
          response = c(response, T)
        } else {
          response = c(response, F)
        }
      }
    }
    names(response) = slot
  }
  response
}

