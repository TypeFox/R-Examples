# Landscape-Methods.R 
# Part of the SEHmodel package.
#
# Copyright (C) 2015        Melen Leclerc <melen.leclerc@rennes.inra.fr>
#                           Jean-Francois Rey <jean-francois.rey@paca.inra.fr>
#                           Samuel Soubeyrand <Samuel.Soubeyrand@avignon.inra.fr>
#                           Emily Walker <emily.walker@avignon.inra.fr>
#                           INRA - BioSP Site Agroparc - 84914 Avignon Cedex 9
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

#' @title simulateInitialPartition Method
#' @description This function creates an object \code{Landscape} and simulates a landscape with neutral and source fields.
#' @name simulateInitialPartition
#' @rdname Landscape-constructor-class
#' @param n Numeric, numbers of fields
#' @param prop Numeric [0,1] toxic fields proportion
#' @param range aggregation parameter (range of the spatial Exponential covariance of Gaussian process)
#' @param xmin x-axis left coordinates
#' @param xmax x-axis right coordinates
#' @param ymin y-axis bottom coordinates
#' @param ymax y-axis top coordinates
#' @return An S4 \code{Landscape} object with \code{n} fields, \code{prop} pourcentage of toxic fields of size (xmin,xmax) (ymin,ymax)
#' @examples
#' # simulate a 5000m x 5000m landscape with 500 fields whose 40% (200) are sources
#'   simulateInitialPartition(n=500,prop=0.4,range=10,xmin=0,xmax=5000,ymin=0,ymax=5000)
#' @export
simulateInitialPartition <- function(n=500,prop=0.4,range=10,xmin=0,xmax=5000,ymin=0,ymax=5000) {
  
  points <- create.voronoi.points(n,prop,range,xmin,xmax,ymin,ymax)
  map <- create.voronoi.diagram(points,xmin,xmax,ymin,ymax)
  
  res <- new("Landscape")
  res@thelandscape <- map
  #proj4string(res@thelandscape)<-CRS("+init=epsg:2154")
  res@xmin <- xmin
  res@xmax <- xmax
  res@ymin <- ymin
  res@ymax <- ymax
  res@n <- n
  
  return(res)
}


#' @title SimulateLandscape Method
#' @name simulateLandscape
#' @description Create an object \code{Landscape}. Simulate a landscape with neutral and source fields and receptors margins.
#' @details Execute both \link{simulateInitialPartition} and \link{simulateThickMargins} functions.
#' @rdname simulateLandscape-constructor-class
#' @param n Numeric, numbers of fields
#' @param prop Numeric [0,1] toxic fields proportion
#' @param range aggregation parameter (range of the spatial exponential covariance of gaussian process) in meters.
#' @param xmin x-axis left coordinates
#' @param xmax x-axis right coordinates
#' @param ymin y-axis bottom coordinates
#' @param ymax y-axis top coordinates
#' @param border_size A numeric, bbox margin
#' @param prob Probability to inflate a filed margin
#' @param mean_thickness Margin width expectation
#' @param v_thickness Margin width variance
#' @return An S4 \code{Landscape} object with n fields, prop pourcentage of toxic fields.
#' @examples \dontrun{
#' simulateLandscape(n=100, prop=0.4, range=10,
#' xmin=0, ymin=1000, xmax=0, ymin=1000, border_size=100,
#' prop=runif(1,0.1,0.9), mean_thickness=runif(1,2,20),
#' v_thickness=50) }
#' @seealso \link{simulateInitialPartition} and \link{simulateThickMargins}
#' @export
simulateLandscape <- function(n=500,prop=0.4,range=10,xmin=0,xmax=5000,ymin=0,ymax=5000, border_size=200,prob=runif(1,0.1,0.9),mean_thickness=runif(1,2,20),v_thickness=50) {
  objectLTemp<-simulateInitialPartition(n=n,prop=prop,range=range,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
  objectL<-simulateThickMargins(objectLTemp,border_size,prob,mean_thickness,v_thickness)
  return(objectL)
}

#' @title Plot Method for Landscape
#' @name plot Landscape
#' @description Plot a Landscape object
#' @param x A Landscape object
#' @param y ANY
#' @param add boolean to add draw hover plot
#' @param time Time selection (default = -1) for Toxic intensity
#' @param objectT Toxic intensity matrix (result of \link{toxicIntensity}, default NULL)
#' @param ... default plot par
#' @param plot.legend plot legend (default TRUE)
# @docType methods
#' @rdname Landscape-plot-method
#' @aliases plot,Landscape,ANY-method
#' @export
setMethod(f="plot",
          signature="Landscape",
          definition=function(x,y,add=F,objectT=NULL,time=-1,...,plot.legend=TRUE){
            alpha=255
            if(add == T) {
              alpha = alpha/3
            }
            else{
              if(plot.legend == TRUE) {
                par(mar=c(5,0,1,5))
                par(xpd=T)
              }
            }
            
            plot(x@thelandscape,add=add,...)
            if(sum(x@thelandscape@data$neutral==1))plot(x@thelandscape[which(x@thelandscape@data$neutral==1),],col=rgb(255,255,255,alpha=alpha,maxColorValue = 255),add=TRUE,...)
            if(sum(x@thelandscape@data$sources==1)) {
              plot(x@thelandscape[which(x@thelandscape@data$sources==1),],col=rgb(255,255,255,alpha=alpha,maxColorValue = 255),add=TRUE,...)
              plot(x@thelandscape[which(x@thelandscape@data$sources==1),],add=TRUE,density=c(9,9),angle=c(-45,45),...)
              plot(x@thelandscape[which(x@thelandscape@data$sources==1),],add=TRUE,density=c(9,9),angle=c(45,-45),...)
            }
            if( "receptors" %in% colnames(x@thelandscape@data) && sum(x@thelandscape@data$receptors==1)) {
              plot(x@thelandscape[which(x@thelandscape@data$receptors==1),],col=rgb(138,43,226,alpha=alpha,maxColorValue = 255),add=TRUE,border=rgb(138,43,226,alpha=alpha,maxColorValue = 255),...)
            }
            if(plot.legend == TRUE){ legend(x=x@thelandscape@bbox[1,1],y=x@thelandscape@bbox[2,1],c("source","neutral","receptor"),pch=c(12,22,22),pt.cex=1.5,pt.lwd=c(2,2,2),col=c(rgb(0,0,0,alpha=alpha,maxColorValue = 255),rgb(0,0,0,alpha=alpha,maxColorValue = 255),rgb(0,0,0,alpha=alpha,maxColorValue = 255)),bg=rgb(255,255,255,alpha=alpha,maxColorValue=255),pt.bg=c(rgb(255,255,255,alpha=alpha,maxColorValue = 255),rgb(255,255,255,alpha=alpha,maxColorValue = 255),rgb(138,43,226,alpha=alpha,maxColorValue = 255)),title="Landscape") }
            
            if(!is.null(objectT) & class(objectT) == "ToxicIntensityRaster" & length(time) == 1 & time > 0) {
              p<-heat.colors(100, alpha = 0.6)
              p[100]=rgb(0,0,0,alpha=0)
              temp<-objectT[time,,]
              temp[which(temp<=0)] <- NA
              r<-raster(as.matrix(temp),crs=CRS("+init=epsg:2154"))
              extent(r)<-extent(x@xmin,x@xmax,x@ymin,x@ymax)
              if( !is.na(proj4string(x@thelandscape)) ) { r<-projectRaster(r,crs=proj4string(x@thelandscape)) }
              raster::image(r,col=p[length(p):1],useRaster=F,add=T,bg="transparent")
              
              if(plot.legend == TRUE) {
                fields::image.plot(as.matrix(temp),legend.only=T,smallplot=c(0.85,0.88,0.20,0.90),col=p[length(p):1])
                mtext(text ="Toxic intensity",line = 0,side = 3,adj = 1.1,padj = 1)
              }
            }
          }
)

# setGeneric(name="simulateLandscape",
#            def=function(object)
#              standardGeneric("simulateLandscape")
# )
# 
# setMethod(f="simulateLandscape",
#           signature="Landscape",
#           definition=function(object){
#             
#           }
# )


#' Show Method for Landscape
#' @param object A Landscape class
#' @rdname Landscape-show-method
#' @aliases show,Landscape-method 
setMethod(f="show",
          signature="Landscape",
          definition=function(object) {
            cat("*** Class Landscape, method Show ***\n")
            cat(sprintf("* @thelandscape : %s\n",class(object@thelandscape)))
            cat(sprintf("* Fields Number = %i\n",object@n))
            if( "sources" %in% colnames(object@thelandscape@data) ) {
              cat(sprintf("  * Sources : %i\n",sum(object@thelandscape@data$sources==1)))
            }
            if( "neutral" %in% colnames(object@thelandscape@data) ) {
              cat(sprintf("  * Neutral : %i\n",sum(object@thelandscape@data$neutral==1)))
            }
            if( "receptors" %in% colnames(object@thelandscape@data) ) {
              cat(sprintf("  * Receptors : %i\n",sum(object@thelandscape@data$receptors==1)))
            }
            cat(sprintf("* xmin = %f\n",object@xmin))
            cat(sprintf("* xmax = %f\n",object@xmax))
            cat(sprintf("* ymin = %f\n",object@ymin))
            cat(sprintf("* ymax = %f\n",object@ymax))
            cat("*** End Show(Landscape) ***\n")  
          }
)

#' Print Method for Landscape
#' @param x A Landscape object
#' @param ... further arguments passed to or from other methods.
#' @rdname Landscape-print-method
#' @aliases print,Landscape-method 
setMethod(f="print",
          signature="Landscape",
          function(x,...) {
            cat("*** Class Landscape, method Print ***\n")
            cat(sprintf("* @thelandscape : %s\n",class(x@thelandscape)))
            cat(sprintf("* Fields Number = %i\n",x@n))
            if( "sources" %in% colnames(x@thelandscape@data) ) {
              cat(sprintf("  * Sources : %i\n",sum(x@thelandscape@data$sources==1)))
            }
            if( "neutral" %in% colnames(x@thelandscape@data) ) {
              cat(sprintf("  * Neutral : %i\n",sum(x@thelandscape@data$neutral==1)))
            }
            if( "receptors" %in% colnames(x@thelandscape@data) ) {
              cat(sprintf("  * Receptors : %i\n",sum(x@thelandscape@data$receptors==1)))
            }
            cat(sprintf("* xmin = %f\n",x@xmin))
            cat(sprintf("* xmax = %f\n",x@xmax))
            cat(sprintf("* ymin = %f\n",x@ymin))
            cat(sprintf("* ymax = %f\n",x@ymax))
            cat("*** End Print(Landscape) ***\n")  
          }
)


#' @title simulateThickMargins Method
#' @name simulateThickMargins
#' @param objectL A Landscape object
#' @param ... other parameters
#' @rdname Landscape-simulateThickMargins-method
#' @exportMethod simulateThickMargins
setGeneric(name="simulateThickMargins",
           def=function(objectL,...)
             standardGeneric("simulateThickMargins")
)

# simulateThickMargins
#' @description Simulate thick margins as receptors in a landscape.
#' @details Margin width use a Gamma distribution with shape and scale parameters based on thickness mean and variance.   
#' @param border_size A numeric, bbox margin
#' @param prob Probability to inflate a margin
#' @param mean_thickness Margin width expectation in meter
#' @param v_thickness Margin width variance in meter
#' @return a Landscape object with fields margin as receptor
#' @aliases simulateThickMargins,Landscape-method
#' @rdname Landscape-simulateThickMargins-method
#' @seealso \link{simulateInitialPartition} and \link{simulateLandscape}
setMethod(f="simulateThickMargins",
          signature="Landscape",
          definition=function(objectL,border_size=200,prob=runif(1,0.1,0.9),mean_thickness=runif(1,2,20),v_thickness=50){
            
            coor<-data.frame("x"=c(objectL@xmin+border_size,objectL@xmax-border_size,objectL@xmax-border_size,objectL@xmin+border_size,objectL@xmin+border_size),
                             "y"=c(objectL@ymax-border_size,objectL@ymax-border_size,objectL@ymin+border_size,objectL@ymin+border_size,objectL@ymax-border_size))
            
            c<-data.frame("x"=c(objectL@xmin-200,objectL@xmax+200,objectL@xmax+200,objectL@xmin-200,objectL@xmin-200),
                          "y"=c(objectL@ymax+200,objectL@ymax+200,objectL@ymin-200,objectL@ymin-200,objectL@ymax+200))
                          
            # set a box
            hole<-Polygon(coor,hole=TRUE)
            o<-Polygon(c,hole=FALSE)
            geom<-Polygons(list(o,hole),"box")
            box<-SpatialPolygons(list(geom))
            
            # Get border
            segments<-edgeslines(objectL@thelandscape) #on récupère les segments (bordures) à partir du pavage
            nsegment<-length(segments@lines)
            
            bernoulli<-rbinom(nsegment,1,prob) # on tire dans la binomiale la proba d'avoir une largeur n
            margin_width<-rgamma(nsegment,shape=mean_thickness*mean_thickness/v_thickness,rate=mean_thickness/v_thickness) #largeur des bordures tirée dans une gamma
            
            length_margin<-bernoulli*margin_width #largeur des bordures de champ
            #length_margin[which(length_margin==0)]<-1e-1
            
            # get polygons ID
            ids=sapply(slot(objectL@thelandscape, "polygons"), function(x) slot(x, "ID"))
            
            #create margins
            seg_margin<-gBuffer(segments,byid=TRUE,width=length_margin/2,capStyle="ROUND") # inflated segments
            
            margin<-gUnionCascaded(seg_margin) # inflated segments union in 1 polygon
            margin@bbox<-objectL@thelandscape@bbox
            field_margin<-gDifference(margin,box,byid=T)
            field_margin@bbox<-objectL@thelandscape@bbox
            proj4string(field_margin)<-proj4string(objectL@thelandscape) 
            margin_ids=as.character(seq(max(as.numeric(ids))+1,max(as.numeric(ids))+length(field_margin@polygons)))
            nb_margins=length(field_margin@polygons)
            
            receptors<-rep(0,nrow(objectL@thelandscape@data)+nb_margins)
            #receptors[(as.numeric(margin_ids)+1)]=1
            receptors[(nrow(objectL@thelandscape@data)+1):(nrow(objectL@thelandscape@data)+nb_margins)]<-1
            
            fields<-gDifference(objectL@thelandscape,field_margin,byid=TRUE,id=ids)
            fields@bbox<-objectL@thelandscape@bbox
            
            # change margins ID
            for(mp in 1:nb_margins) {
              field_margin@polygons[[mp]]@ID=margin_ids[mp]
            }
            # Add margin to polygons landscape
            temp_landscape<-SpatialPolygons(c(fields@polygons,field_margin@polygons))
            
            # Add data.frame with receptors added
            new_data_frame<-cbind(rbind(objectL@thelandscape@data,data.frame("sources"=rep(0,nb_margins),"neutral"=rep(0,nb_margins),row.names = margin_ids[])),"receptors"=receptors)
            new_landscape<-SpatialPolygonsDataFrame(temp_landscape, new_data_frame)
            new_landscape@bbox<-objectL@thelandscape@bbox
            proj4string(new_landscape)<-proj4string(objectL@thelandscape)
            
            objectL@thelandscape<-new_landscape
            
            return(objectL)
          }
)

# Method getSpatialPolygons
# @name getSpatialPolygons
# @rdname Landscape-getSpatialPolygons-method
# @exportMethod getSpatialPolygons
setGeneric(name="getSpatialPolygons",
           def=function(object,...)
             standardGeneric("getSpatialPolygons")
)

# Get all polygons of a landscape object
# @param object A Landscape object
# @return a SpatialPolygonDataFrame object
setMethod(f="getSpatialPolygons",
          signature="Landscape",
          definition=function(object) {
            return(object@thelandscape)            
          }
)

#' Method to get SpatialPolygons Sources from a Landscape
#' @name getSPSources
#' @param object A Landscape object
#' @param ... options
#' @examples 
#' \dontrun{
#' ## To get the first spatial polygon source information from a landscape object (named land)
#' getSPSources(land)@polygons[[1]]
#' }
#' @rdname Landscape-getSPSources-method
#' @exportMethod getSPSources
setGeneric(name="getSPSources",
           def=function(object,...)
             standardGeneric("getSPSources")
)
 
# Method to get SpatialPolygons Sources from a Landscape
# @name GetSPSources
#' @description Get all polygons of a Landscape object identified as sources
#' @aliases Landscape,getSPSources-method
#' @return a SpatialPolygonDataFrame object
#' @rdname Landscape-getSPSources-method
setMethod(f="getSPSources",
          signature="Landscape",
          definition=function(object) {
            return(object@thelandscape[which(object@thelandscape@data$sources==1),])            
          }
)

# Method getSPReceptors
# @name getSPReceptors
# @rdname Landscape-getSPReceptors-method
# @exportMethod getSPReceptors
setGeneric(name="getSPReceptors",
           def=function(object,...)
             standardGeneric("getSPReceptors")
)
# @name getSPReceptors
#
# Get all polygons of a landscape object identify as receptors
#
# @param object A Landscape object
# @return a SpatialPolygonDataFrame object
# @aliases Landscape,getSPReceptors-method
setMethod(f="getSPReceptors",
          signature="Landscape",
          definition=function(object) {
            return(object@thelandscape[which(object@thelandscape@data$receptors==1),])            
          }
)

#' Wrapper function : loadLandscape
#' @name loadLandscape
#' 
#' @description Wrapper function to create a Landscape object using SpatialPolygons and dataframe.
#' The SpatialPolygons object and the data.frame have to contain the same number of polygons
#'  and row (row ID is polygons ID).
#' 
#' @rdname Landscape-load-class
#' @param sp a SpatialPolygons object designing the landscape
#' @param data a data.frame containing fields (polygons) information. Row names as fields ID, column names as sources | neutral | receptors (for a given field, the value is 1 for the type of the field (source or neutral or receptor), otherwise 0).
#' @return A Landscape object
#' @export
loadLandscape <- function(sp,data) {
  res <- new("Landscape")
  
  sptemp<-spTransform(sp,CRSobj = CRS("+init=epsg:2154"))
  newsp<-SpatialPolygons(sptemp@polygons)
  proj4string(newsp)<-proj4string(sptemp)
  res@thelandscape <- SpatialPolygonsDataFrame(newsp,data = data)
  res@thelandscape@bbox<-sptemp@bbox
  res@thelandscape@plotOrder<-sptemp@plotOrder
  res@xmin <- extent(sptemp)[1]
  res@xmax <- extent(sptemp)[2]
  res@ymin <- extent(sptemp)[3]
  res@ymax <- extent(sptemp)[4]
  res@n <- length(sp@polygons)
  
  return(res)
  
}

#' Create a Landscape object from SIG shapefile file
#' @name loadLandscapeSIG
#' 
#' @description Create a Landscape object from SIG shapefile. Shapefile has to contain a SpatialPolygonsDataFrame. Data in the data frame contain fields (polygons) information. Row names as fields ID, cols names as sources | neutral | receptors (for a given field, the value is 1 for the type of the field (source or neutral or receptor), otherwise 0).
#' 
#' @rdname Landscape-load-sig-class
#' @param dsn folder path to the source files
#' @param layer file name without extension
#' @return A Landscape object
#' @export
loadLandscapeSIG <- function(dsn,layer) {
  
  if(requireNamespace("rgdal",quietly=T)) {
  ogr <- readOGR(dsn,layer)
  res <- new("Landscape")
  
  res@thelandscape<-spTransform(ogr,CRSobj = CRS("+init=epsg:2154"))
  res@n <- length(res@thelandscape@polygons)
  res@xmin <- res@thelandscape@bbox[1]
  res@xmax <- res@thelandscape@bbox[3]
  res@ymin <- res@thelandscape@bbox[2]
  res@ymax <- res@thelandscape@bbox[4]
  return(res)
  }
  else {
    warning("loadLandscapeSIG function need \"rgdal\" package")
    return(NULL)
  }
  
}


#' Save Particles Dispersion 3D Array to tiff file
#' @name saveIntoFile
#' 
#' @description Save into tiff file particles dispersion 3D array from toxicIntensity.
#' The output is a RasterStack with a layer per time unit with projection set to CRS="+proj=longlat +datum=WGS84"
#' @rdname Landscape-save-tiff
#' @param objectL a Landscape object
#' @param objectT a 3D array particles dispersion indexed by time
#' @param filename output file name
#' @param format output format (default=GTiff)
#' @param overwrite if True overwrite filename
#' @return a RasterStack object
#' @export
saveIntoFile<-function(objectL,objectT,filename="ParticlesDispersion.tif",format="GTiff",overwrite=T) {
  
  timemax=length(objectT[,1,1])
  
  r<-raster(as.matrix(objectT[1,,]),crs=objectL@thelandscape@proj4string)
  #extent(r)<-extent(objectL@thelandscape)
  names(r)<-paste("time",as.character(1),sep=".")
  
  for(t in 2:timemax) {
    rtemp<-raster(as.matrix(objectT[t,,]),crs=objectL@thelandscape@proj4string)
    #extent(rtemp)<-extent(objectL@thelandscape)
    names(rtemp)<-paste("time",as.character(t),sep=".")
    r<-addLayer(r,rtemp)
  }
  extent(r)<-extent(objectL@thelandscape)
  
  r<-projectRaster(r,crs="+proj=longlat +datum=WGS84")
  
  rf <- writeRaster(r, filename=filename, format=format, overwrite=overwrite)  #bylayer=TRUE
  
  return(rf)
  
}

#################
#### PRIVATE ####
#################

### create.voronoi.points : fonction qui simule les graines pour le pavage de voronoi avec des marques OGM NonOGM
## pour contrôler l'agregation spatiale des points OGM on utilise un processus gaussien avec autocorrélation spatiale
## ici on utilise une covariance exponentielle. 
# L'idée étant de faire une analyse de sensibilité on fixe le paramètre phi et on fait varier l'étendue r
##prend en paramètre le nombre de graines n, la proportion ogm prop, xmin xmax ymin ymax, le paramètre de porté de la covariance covar
create.voronoi.points <- function(n,prop,range,xmin,xmax,ymin,ymax){
  
  coor<-cbind(runif(n,xmin,xmax),runif(n,ymin,ymax))
  d<-as.matrix(dist(coor))
  cov_mat<-Exponential(d,range=range,phi=10)
  s<-mvrnorm(1,mu=rep(0,n),Sigma=cov_mat)
  n_ogm=floor(prop*n) #nombre de parcelles OGM
  data<-as.data.frame(s)
  sources<-rep(0,n)
  neutral<-rep(0,n)
  threshold<-max(sort(data$s)[1:n_ogm])
  sources[which(s<=threshold)]=c(1) #marque source (OGM) = 1 dans la colonne sources du dataframe
  neutral[which(s>threshold)]=c(1) #marque source neutre = 1 dans la colonne neutral du dataframe
  f<-SpatialPointsDataFrame(coords=coor,data=cbind.data.frame(sources,neutral))
  return(f)
}


## vornoi retourne un objet SpatialPolygonsDataFrame correspondant à un pavage de voronoi
#layer= SpatialPointsDataFrame, xmin,xmax, ymin, ymax= bornes du domaine/paysage
create.voronoi.diagram <- function(layer,xmin,xmax,ymin,ymax){
  crds = layer@coords     #extraction de la matrice des coordonn?es de l'objet layer
  z = deldir(crds[,1], crds[,2],rw=c(xmin,xmax,ymin,ymax),suppressMsge=T) # triangulation en fonction des points proposés
  w = tile.list(z)      #pour chaque point-centre extrait le polygone-pavage
  polys = vector(mode='list', length=length(w))
  for (i in seq(along=polys)) {
    pcrds = cbind(w[[i]]$x, w[[i]]$y)
    pcrds = rbind(pcrds, pcrds[1,])
    polys[[i]] = Polygons(list(Polygon(pcrds)), ID=as.character(i))
  }
  SP = SpatialPolygons(polys)
  SP@bbox[,1]<-c(xmin,ymin)
  SP@bbox[,2]<-c(xmax,ymax)
  voronoi = SpatialPolygonsDataFrame(SP,data=layer@data)
  return (voronoi)
}

## edgeslines : à partir d'un objet polygons, retourne un SpatialLines (mieux ordonné qu'avec as(,SpatialLines))
edgeslines = function(voronoi)
{
  mat = matrix(c(0,0,0,0),ncol=4,byrow=TRUE)
  for (i in seq(voronoi))
  {
    coor = voronoi@polygons[[i]]@Polygons[[1]]@coords
    for (j in seq(2,nrow(coor)))
    {
      if (coor[j-1,1]<coor[j,1]) { mat = rbind(mat,c(coor[j-1,],coor[j,])) }
      else { mat = rbind(mat,c(coor[j,],coor[j-1,])) }
    }    
  }
  mat = unique(mat[-1,])
  edges = vector(mode='list',length=nrow(mat))
  for (i in seq(nrow(mat))) { edges[[i]] = Lines(Line(coords=rbind(mat[i,1:2],mat[i,3:4])),ID=as.character(i)) }
  SE = SpatialLines(edges)
}

