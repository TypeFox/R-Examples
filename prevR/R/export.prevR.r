#' @exportMethod export
setGeneric("export",
    function(object, element, format, file ,N= NULL, R = NULL, clusters.only = FALSE, ext=NULL, sep=NULL, dec=NULL,  ...){ 
        standardGeneric("export") 
    }
)

#' Export an object of class prevR.
#' 
#' This method could be used to export an object of class \code{\link[=prevR-class]{prevR}} 
#' in different formats (text, shapefile, dbase...)
#' 
#' @param object object of class \code{\link[=prevR-class]{prevR}}.
#' @param element element to export: "clusters" or "boundary".
#' @param format format: "dbf", "txt", csv", "csv2" or "shp" (unused if \code{element="boundary"}).
#' @param file file name without extension.
#' @param N integer or list of integers setting elements of \code{rings} to export 
#'   (unused if \code{element="boundary"}).
#' @param R integer or list of integers setting elements of \code{rings} to export 
#'   (unused if \code{element="boundary"}).
#' @param clusters.only export only the slot \code{clusters} of \code{object} 
#'   (unused if \code{element="boundary"})?
#' @param ext coerce the extension of the export file 
#'   (unused if \code{element="boundary"} or if \code{format="shp"}).
#' @param sep coerce the field separator string (unused if \code{element="boundary"} 
#'   or if \code{format="shp"} or if \code{format="dbf"}).
#' @param dec coerce the string to use for decimal point 
#'   (unused if \code{element="boundary"} or if \code{format="shp"} or if \code{format="dbf"}).
#' @param \dots additional arguments transmitted to \code{\link[maptools]{writePolyShape}}, 
#'   \code{\link[maptools]{writePointsShape}}, \code{\link[foreign]{write.dbf}} or 
#'   \code{\link[utils]{write.table}}.
#' 
#' @details If \code{element="boundary"}, the slot \code{boundary} of \code{object} 
#' will be exported as a \emph{shapefile}.
#' 
#' Otherwise, the slot \code{clusters}, merged with the slot \code{rings}, will be exporter.\cr
#' See \code{\link{as.data.frame.prevR}} for details on the use of the parameters of \code{N}, 
#' \code{R} et \code{clusters.only}.
#' 
#' \code{format} specifies the export format of the data frame returned by 
#' \code{\link{as.data.frame.prevR}}: \tabular{ll}{
#'     "shp" \tab Shape File (require packages \pkg{maptools} and \pkg{foreign})\cr
#'     "dbf" \tab DBASE format (extension: .dbf, require the package \pkg{foreign})\cr
#'     "txt" \tab tabulated text (extension: .txt)\cr
#'     "csv" \tab 'comma separated values' (extension: .csv)\cr
#'     "csv2" \tab CSV variant using a semicolon as field separator (extension: .csv)
#' }
#' \code{ext} could be used to coerce the extension of the output filen except for 
#' \emph{shapefile} export, which will write three different files (.shp, .shx et .dbf).
#' 
#' The "txt" format uses by default a tabulation as field separator and a point "." for decimal point. 
#' The "csv" format uses a comma "," as field separator and a point "." as decimal point.\cr
#' The "csv2" format is a variant using a semicolon ";" as field separator and a colon "," for decimal point, 
#' the Excel convention for CSV files in some Western European locales. \cr
#' \code{sep} and \code{dec} could be used to coerce the field separator and the decimal point 
#' (together with the "txt" format).
#' 
#' @seealso \code{\link[maptools]{writePolyShape}} \{\pkg{maptools}\}, \code{\link[maptools]{writePointsShape}}
#' \{\pkg{maptools}\}, \code{\link[foreign]{write.dbf}} \{\pkg{foreign}\}, \code{\link[utils]{write.table}}
#' \{\pkg{utils}\}.
#' 
#' @examples 
#'   \dontrun{
#'     export(fdhs, element="boundary", file="area")
#'     export(fdhs, element="clusters", format="shp", file="points")
#'     
#'     dhs <- rings(fdhs,N=c(100,300,500))
#'     export(dhs, element="clusters", format="csv", N=300, file="points")
#'   }
#' 
#' @aliases export-methods export,prevR-method export
#' @keywords manip spatial

setMethod("export","prevR",
  function(object, element, format, file, N= NULL, R = NULL, clusters.only = FALSE, ext=NULL, sep=NULL, dec=NULL, ...){
  ##################################################################################################
  # Cette fonction exporte dans des fichiers de differents formats les elements d'un objet de class prevR
  # les differents arguments sont
  #    object : un objet de classe prevR
  #    element   : Une chaine character contenant "clusters" ou "boundary"
  #           Si element = clusters on exporte un merge de clusters et des rings definis par N et R 
  #                    (Dans le cas ou l'argument clusters.only est positionne a T on n'exporte que clusters)
  #           Si element = boundary on exporte au format shp l'element boundary
  #   format  : une chaine de character specifiant le format du fichier exporte
  #            Cet argument n'est actif que si element = clusters
  #            Si format = dbf on exporte au format dbf
  #            Si format = txt on exporte en texte tabule (point pour les decimales)
  #            Si format = csv on exporte en csv (separateur virgule et point pour les decimales)
  #            Si format = csv2 on exporte en csv (separateur point-virgule et virgule pour les decimales)
  #            Si format = shp on exporte au format shp
  #   file  : le nom du fichier sans son extension
  #           Si le format est shp 3 fichiers aux extensions  *.shp, *.shx et *.dbf   sont crees
  #   ext :    permet de forcer l'extension du fichier cree (sans effet pour les formats shapefiles)
  #               si NULL, extension dbf pour le format dbf et txt pour le format txt.
  #   ...  : permettent de passer tous les arguments des fonctions d'ecritures
  #    
  # Les fonctions d'ecriture appelees sont
  #      writePolyShape pour l'ecriture de l'element boundary (package maptools)
  #      writePointsShape pour l'ecriture de l'element cluster (package maptools)
  #      write.table
  #      write.dbf (package foreign)
  # 
  # 
  # Exemples
  #    export(data.prevR, element="clusters", format="text", file="c:/Temp/BF6.csv",sep = ";",  clusters.only = F)
  #    export(data.prevR, element="clusters", format="shp", file="c:/Temp/BF7",clusters.only = F)
  #    export(data.prevR, element="boundary", file="c:/Temp/BF8")     
  ##################################################################################################
  # Dans ... on passe les arguments de as.data.frame
    ind = match(element,c("clusters","boundary"),nomatch=0)
    if(ind == 0){
      stop("the 'element' argument must be 'clusters' or 'boundary'.", call.=F) 
    }
    if(element=="boundary"){
      boundary = slot(object,"boundary")
      if(attr(boundary,"valid")){
        IDs <- sapply(slot(boundary, "polygons"), function(x) slot(x, "ID"))
        data <- data.frame(1, IDs,row.names=IDs)
        names(data) <- c("id", "name")
        SPDF <- SpatialPolygonsDataFrame(boundary, data)
        maptools::writePolyShape(SPDF, file, ...)
      }
      return(NULL)
    }
    
    if(element=="clusters"){
      ind = match(format,c("dbf","txt","shp","csv","csv2"),nomatch=0)
      if(ind == 0){
        stop("the 'format' argument must be 'dbf', 'txt', 'csv', 'csv2' or 'shp'.", call.=F) 
      }
      clusters = as.data.frame(object, N = N , R = R, clusters.only = clusters.only)
      if(format=="shp") {
        PS = clusters
        coordinates(PS)= ~x+y
        PS@proj4string = object@proj
        maptools::writePointsShape(PS,file, ...)
      }
      if (is.null(ext) && format!='csv2') ext = format
      if (is.null(ext) && format=='csv2') ext = 'csv'
      file = paste(file,ext,sep='.')
      if(format=="txt") {
        if (is.null(sep)) sep="\t"
        if (is.null(dec)) dec="."
        write.table(clusters, file = file,row.names = F,sep = sep, dec = dec, ...)
      }
      if(format=="csv") {
        if (is.null(sep)) sep=","
        if (is.null(dec)) dec="."
        write.table(clusters, file = file,row.names = F,sep = sep, dec = dec, ...)
      }
      if(format=="csv2") {
        if (is.null(sep)) sep=";"
        if (is.null(dec)) dec=","
        write.table(clusters, file = file,row.names = F,sep = sep, dec = dec, ...)
      }
      if(format=="dbf")  {
        foreign::write.dbf(clusters, file, ...)
      }
    }
  }
)

       