#' @exportMethod changeproj

setGeneric("changeproj",
    function(object,proj){ 
        standardGeneric("changeproj") 
    }
)

#' Convert map projection of a object of class prevR.
#' 
#' This function converts map projection (and/or datum) used by an object of 
#' class \code{\link[=prevR-class]{prevR}} into another one.
#' 
#' @param object object of class \code{\link[=prevR-class]{prevR}}.
#' @param proj new map projection.
#' 
#' @details \code{proj} could be a character string corresponding to a 
#' \emph{PROJ.4} projection (see \url{http://trac.osgeo.org/proj/} for more details) 
#' or an object of class \code{\link[sp:CRS-class]{CRS}}\{\pkg{sp}\}.
#' 
#' \code{changeproj} transform the columns  "x" and "y" of the slot \code{clusters} of 
#' \code{object} and convert \code{boundary} using the new map projection defined by \code{proj}.\cr
#' If applicable, the slot \code{rings} will be recalculated.
#' 
#' @return Return \code{object} expressed in the projection \code{proj}.
#' @seealso \code{\link[rgdal]{spTransform}}\{\pkg{rgdal}\}, \code{\link{prevR-class}}.
#' 
#' @examples 
#' print(fdhs)
#' plot(fdhs, axes=TRUE, main="Projection: longitude/latitude")
#' 
#' fdhs2 <- changeproj(fdhs,
#'                    "+proj=utm +zone=30 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
#' print(fdhs2)
#' dev.new()
#' plot(fdhs2, axes=TRUE, main="Projection: UTM Zone 30")
#' 
#' @keywords manip spatial
#' @aliases changeproj changeproj-methods
setMethod("changeproj","prevR",
  function(object, proj){
  ###############################################################################################  
  # cette fonction transforme un objet de classe prevR dans une autre projection definie par proj 
  # Cette fonction modifie donc
  #   Les colonnes x et y du slot clusters
  #   Le slot boundary
  # Pour faire cela on utilise la fonction spTransform du package rgdal
  # proj peut etre une chaine de caractere ou un objet CRS
  ###############################################################################################
    if (class(proj)!="CRS")
      proj = CRS(proj)
    
    # cluster slot modification
    clusters                = slot(object,"clusters")
    coordinates(clusters)   = c("x", "y")
    proj4string(clusters)   = slot(object,"proj")
    clusters                = spTransform(clusters, proj)
    xy                      = slot(clusters,"coords")
    clusters                = slot(object,"clusters")
    clusters[,c("x","y")]   = xy
    slot(object,"clusters") = clusters
    
    # boundary slot modification
    boundary                = slot(object,"boundary")
    is.valid                = attr(boundary,"valid")
    proj4string(boundary)   = slot(object,"proj")
    boundary                = spTransform(boundary,proj)
    attr(boundary,"valid")  = is.valid
    
    N = NULL
    R = NULL
    # On verifie si rings est vide
    if(is.prevR(object,"rings")){
      rings = slot(object,"rings")
      N = sapply(rings,function(x) x$N)
      R = sapply(rings,function(x) x$R)
    }
    
    slot(object,"boundary") = boundary
    slot(object,"proj")     = proj
    slot(object,"rings")    = list()
    
    # On recalcule, le cas echeant le slot rings
    if(!is.null(N) && !is.null(R))
      object = rings(object,N=N,R=R)
    
    object
  }
)

