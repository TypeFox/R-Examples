#' Objects of class prevR.
#' 
#' Class used by the package \pkg{prevR}
#' 
#' @section Objects from the Class:
#' Objects of this class could be created by the function \code{\link{as.prevR}}.
#' @slot clusters \code{data.frame} with observed data (one line per cluster). 
#' Columns names are:\itemize{
#'   \item "id" cluster ID.
#'   \item "x" longitude.
#'   \item "y" latitude.
#'   \item "n" number of valid observations per cluster.
#'   \item "pos" number of positive cases per cluster.
#'   \item "prev" observed prevalence (in \%) in the cluster (pos/n).
#'   \item "wn" (optional) sum of weights of observations per cluster.
#'   \item "wpos" (optional) sum of weights of positive cases per cluster.
#'   \item "wprev" (optional) weighted observed prevalence (in \%) in the cluster (wpos/wn).
#'   \item "c.type" (optional) cluster type.
#' }
#' @slot boundary object of class \code{\link[sp:SpatialPolygons-class]{SpatialPolygons}}, 
#' borders of the studied area.
#' @slot proj object of class \code{\link[sp:CRS-class]{CRS}}, map projection used.
#' @slot rings list of results returned by \code{\link[=rings,prevR-method]{rings}}. 
#' Each entry is composed of 3 elements: \code{N}, minimum number of observations per ring; 
#' \code{R}, maximum radius of rings and \code{estimates}, a data frame with the following 
#' variables:\itemize{
#'   \item "id" cluster ID.
#'   \item "r.pos" number of positive cases inside the ring.
#'   \item "r.n" number of valid observations inside the ring.
#'   \item "r.prev" observed prevalence (in \%) inside the ring (r.pos/r.n).
#'   \item "r.radius" ring radius (in kilometers if coordinates in decimal degrees, 
#'     in the unit of the projection otherwise).
#'   \item "r.clusters" number of clusters located inside the ring.
#'   \item "r.wpos" (optional) sum of weights of positive cases inside the ring.
#'   \item "r.wn" (optional) sum of weights of valid observations inside the ring.
#'   \item "r.wprev" (optional) weighted observed prevalence (in \%) inside the ring (r.wpos/r.wn).
#' }
#' Note: the list \code{rings} is named, the name of each element is 
#' N\emph{N_value}.R\emph{R_value}, for example \emph{N300.RInf}.
#' @section Methods: \describe{
#'   \item{as.data.frame}{\code{signature(x = "prevR")} converts an object of class prevR 
#'     into a data frame.}
#'   \item{as.SpatialGrid}{\code{signature(object = "prevR")} generates a spatial grid.}
#'   \item{export}{\code{signature(object = "prevR")} exports a prevR object as a shapefile, 
#'     a dbase file or a text file.}
#'   \item{idw}{\code{signature(formula = "ANY", locations = "prevR")} calculates a spatial 
#'     interpolation using an inverse distance weighting.}
#'   \item{kde}{\code{signature(object = "prevR")} estimates a prevalence surface using 
#'     kernel density estimators.}
#'   \item{krige}{\code{signature(formula = "ANY", locations = "prevR")} calculates a 
#'     spatial interpolation by kriging.}
#'   \item{plot}{\code{signature(x = "prevR", y = "ANY")} plots data of a prevR object.}
#'   \item{print}{\code{signature(x = "prevR")} shows a summary of a prevR object.}
#'   \item{rings}{\code{signature(object = "prevR")} calculates rings of equal number of 
#'     observations and/or equal radius.}
#'   \item{show}{\code{signature(object = "prevR")} shows a summary of a prevR object.}
#'   \item{summary}{\code{signature(object = "prevR")} shows a summary of the variables 
#'     of a prevR object.}
#'   \item{changeproj}{\code{signature(object = "prevR")} changes the map projection used.}
#' }
#' @seealso 
#' \code{\link{as.prevR}}, \code{\link{is.prevR}}, \code{\link{changeproj,prevR-method}}, 
#' \code{\link{rings,prevR-method}}, \code{\link{print,prevR-method}}, 
#' \code{\link{plot,prevR-method}}, \code{\link{summary,prevR-method}}, 
#' \code{\link{kde,prevR-method}}, \code{\link{krige,prevR-method}}, 
#' \code{\link{idw,prevR-method}}, \code{\link{export,prevR-method}}.
#' @examples 
#' showClass("prevR")
#' 
#' col <- c(id = "cluster", 
#'         x = "x",
#'         y="y",
#'         n="n",
#'         pos = "pos",
#'         c.type = "residence",
#'         wn="weighted.n",
#'         wpos="weighted.pos"
#' )
#' dhs <- as.prevR(fdhs.clusters,col, fdhs.boundary)
#' str(dhs)
#' print(dhs)
#' 
#' \dontrun{
#'  dhs <- rings(fdhs,N=c(100,300,500))
#'  str(dhs)
#'  print(dhs)
#' }
#' @keywords classes
#' @export

setClass(Class = "prevR",
  ###############################################################################################
  # Definition de la classe prevR
  # un objet prevR a 4 slots
  # slot clustres
  #    un data frame contenant les colonnes
  #        "id" Identifiant du cluster
  #        "x"  Longitude ou toute projection en x du cluster
  #        "y"  Longitude ou toute projection en y du cluster
  #        "n"  Le nombre de personnes entrant dans l'enquete pour ce cluster
  #        "pos" Le nombre de cas positifs pour ce cluster
  #        "wn" (facultatif) idem n mais pondere 
  #        "wpos" (facultatif) idem pos mais pondere
  #        "c.type" (facultatif) Une variable categorielle fournissant le type de chaque cluster 
  #                 Exemple c.type peut contenir le type de localisation (Urban, not Urban) du cluster
  #                 Cette information n'est utilisee que dans la fonction plot quand l'argument type = "c.type"
  #        "prev" prevalence observee 100*pos/n
  #        "wprev" (facultatif) prevalence observee ponderee 100*wpos/wn
  # slot boundary 
  #    un objet de class spatialPolygons contenant en general les frontieres d'un pays
  # slot proj
  #    une chaine de character contenant la projection dans laquelle sont exprimees les donnees
  #        c'est a dire les colonnes x et y de clusters et les polygones de boundary
  # slot rings  (Ce slot est automatiquement rempli par la fonction rings)
  #    Une liste de rings 
  # Chaque ring est une liste contenant 3 elements
  #    N : contient la valeur de N pour laquelle ce ring a ete calcule
  #    R : contient la valeur de R pour laquelle ce ring a ete calcule
  #    estimates : contient un dataframe dont les colonnes sont
  #       id : l'identificateur du cercle (qui est egal a l'identificater du cluster centre du cercle)
  #       r.pos : le nombre de cas positifs dans le cercle
  #       r.n : l'effectif du cercle
  #       r.prev : la prevalence du cercle (100*r.pos/r.n)
  #       r.radius :  le rayon du cercle
  #       r.clusters : le nombre de clusters dans le cercle
  #       Si des donnees ponderes sont presentes 
  #       r.wpos : le nombre pondere de cas positifs dans le cercle
  #       r.wn : l'effectif pondere du cercle
  #       r.wprev : la prevalence ponderee du cercle (100*r.wpos/r.wn)
  # Remarque la list rings est nommee. Le nom de chacun de ses elements est de la forme  
  #  N(valeur de n ).R(valeur de R) : exemple N500.RInf    
  # 
  # 
  ###############################################################################################      
representation(clusters = "data.frame", boundary = "SpatialPolygons", proj = "CRS", rings = "list"),
validity = function(object){
  clusters = slot(object,"clusters")
  necessaryVar = c("id","x","y","n","pos")
  ind = match(necessaryVar,names(clusters))
  if(any(is.na(ind))) {
    missing.var = paste(necessaryVar[is.na(ind)],collapse=", ")
    n.missing = length(necessaryVar[is.na(ind)])
    sprintf(ngettext(n.missing,"the variable %s is missing.","the following variables (%s) are missing.",domain="R-prevR"),missing.var) -> stop.mess
    stop(stop.mess, call.=F)
  }
  coupledVar = c("wn","wpos")
  ind = match(coupledVar,names(clusters))
  ind = ind[!is.na(ind)]
  if(length(ind) == 1){
    stop(gettextf("the wn and wpos variables are coupled and you have only defined %s.",names(clusters)[ind],domain="R-prevR"), call.=F)
  }
  proj =  slot(object,"proj")
  isOk = try(CRS(proj@projargs),silent=T)
  if(attr(isOk,"class") == "try-error"){
    stop(gettextf("the projection %s, defined in the 'proj' argument, is incorect.",proj,domain="R-prevR"), call.=F)
  }
  boundary = slot(object,"boundary")
  boundary.proj = slot(boundary,"proj4string")
  if (boundary.proj@projargs != proj@projargs) {
    stop("'boundary' and 'clusters' didn't have the same projection.", call.=F)
  }
  T
}
)
