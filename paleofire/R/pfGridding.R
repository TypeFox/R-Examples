## 

#' Produce gridded maps of transformed charcoal values.
#' 
#' The function uses weighted spatio-temporal interpolation to produce gridded
#' maps of transformed charcoal values. Spatial grids are used to interpolate
#' transformed charcoal values for a key period defined by \code{Age}. For each
#' grid cell the function search charcoal sites located in a radius defined by
#' \code{distance_buffer} from the grid centre and at an elevation within a
#' range defined by \code{elevation_buffer} from the mean elevation of the
#' cell. Then the function search for charcoal samples within a temporal range
#' from the key date defined by \code{time_buffer}. Finally a tricube distance
#' weighting function is applied to each sample by considering it spatial
#' distance to the grid centre and it temporal distance to the key date. This
#' approach that weight samples according to their spatio-temporal location
#' also down-weight charcoal sites that are poorly sampled.
#' 
#' 
#' @param data An object returned by \code{\link{pfTransform}}.
#' @param cell_sizex Numeric, grid cell width (m).
#' @param cell_sizey Numeric, grid cell height (m).
#' @param age Numeric, key date (years BP).
#' @param cell_size Numeric, grid cell size (bypass cell_sizex and cell_sizey
#' and produce square cells).
#' @param time_buffer Numeric, temporal distance (years) from the key date to
#' search for charcoal samples.
#' @param distance_buffer Numeric, spatial distance from the grid centres to
#' search for charcoal samples (m).
#' @param raster_extent Numeric, define custom extent for the analysis such as
#' raster_extent = c(xmin, xmax, ymin, ymax)
#' @param elevation_buffer Numeric, elevation range from the mean grid cell
#' elevation to search for charcoal sites.
#' @param proj4 String, proj.4 string representing the desired projection for
#' plotted maps. Default is Robinson ("+proj=robin +lon_0=0 +x_0=0 +y_0=0
#' +ellps=WGS84 +datum=WGS84 +units=m +no_defs"). See
#' \url{http://www.spatialreference.org} to look up the string for your
#' favorite projections.
#' @param sea_mask Logical, mask cells falling in the sea.
#' @param other_mask A sp object (SpatialPolygonsDataFrame) used to mask data
#' i.e. for not interpollating pixels under the mask (classical usage: ice
#' extent mask). Note that the SpatialPolygonsDataFrame projection must be used
#' in the analysis and defined using \code{proj4} argument, otherwise the mask
#' should be reprojected (e.g. using rgdal::spTransform).
#' @param verbose Logical, verbose or not...
#' @return A "pfGridding" object (list) that could be plotted using
#' \code{\link{plot.pfGridding}}.
#' @author O.Blarquez
#' @seealso \code{\link{plot.pfGridding}}, \code{\link{pfTransform}},
#' \code{\link{pfDotMap}}
#' @examples
#' 
#' ID=pfSiteSel(id_region=="ENA0", l12==1, long>-85)
#' 
#' TR=pfTransform(ID,method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,4000))
#' 
#' p=pfGridding(TR,age=1000)
#' summary(p)
#' 
#' \dontrun{
#' require(raster)
#' plot(p$raster)
#' 
#' ## Example of other_mask usage: we will use here Dyke 2003 ice extent map for North 
#'  America
#' require(maptools)
#' ID=pfSiteSel(id_region=="ENA0", long>-100,lat>40)
#' TR=pfTransform(ID,method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,4000))
#' 
#' ## Define projection (same as Dyke 2003)
#' proj4="+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 
#'  +lon_0=-95 +x_0=0 +y_0=0 +ellps=clrk66 +datum=NAD27 +units=m +no_defs"
#' 
#' ## Download the shapefile
#' where=getwd() 
#' download.file("http://blarquez.com/public/data/ice_9500_calBP_lcc.shp",
#'                paste0(where,"/ice_9500_calBP_lcc.shp"))
#' download.file("http://blarquez.com/public/data/ice_9500_calBP_lcc.dbf",
#'                paste0(where,"/ice_9500_calBP_lcc.dbf"))
#' download.file("http://blarquez.com/public/data/ice_9500_calBP_lcc.shx",
#'                paste0(where,"/ice_9500_calBP_lcc.shx"))
#' 
#' ice_shp=readShapePoly(paste0(where,"/ice_9500_calBP_lcc.shp"), 
#'                       proj4string=CRS(proj4))
#' plot(ice_shp) 
#' 
#' p=pfGridding(TR,age=9500,cell_size=100000,distance_buffer=300000,
#'              proj4=proj4,other_mask=ice_shp)
#' plot(p,add=ice_shp)
#' 
#' # Citation: Dyke, A.S., Moore, A. And Robertson, L. 2003 :
#' # Deglaciation of North America, Geological Survey of Canada Open File 1574.
#' }
#' 
pfGridding=function(data,cell_sizex=NULL,
                    cell_sizey=NULL,
                    age=0,
                    cell_size=NULL,
                    time_buffer=NULL,
                    distance_buffer=NULL,
                    raster_extent=NULL,
                    elevation_buffer=NULL,
                    proj4=NULL,
                    sea_mask=TRUE,
                    other_mask=NULL,
                    verbose=TRUE){
  
  ## pfTransform object
  if(class(data)=="pfTransform"){
    data<-data.frame(
      x=rep(summary(data$params$ID)$long,each=length(data$TransData[,1])),
      y=rep(summary(data$params$ID)$lat,each=length(data$TransData[,1])),
      age=c(data$Age),
      char=c(data$TransData));
    data=na.omit(data)
    ## Remove extreme outliers i.e. 3*sd (sometimes happens... why? probably baseperiod related)
    data=data[data[,4]<3*sd(data[,4]) & data[,4]>-(3*sd(data[,4])) ,]
  }
  
  
  ## Limit the input to age+-time_buffer
  if(is.null(time_buffer)) time_buffer=500
  data=data[data[,3]>=age-time_buffer & data[,3]<=age+time_buffer,]
  
  
  # source("/Users/Olivier/Documents/BorealTreeCover/final/triCube.R")
  ## Load countries with lakes from http://www.naturalearthdata.com/downloads/10m-cultural-vectors/
  # load(file="/Users/Olivier/Documents/BorealTreeCover/final/world_map.rda")
  longlat="+proj=longlat|+proj = longlat|+proj =longlat|+proj= longlat|+proj=longlat"
  
  if(is.null(proj4))
    proj4<-"+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  
  xy <-cbind(data[,1],data[,2])
  
  ### LOAD DEM if required
  if(is.null(elevation_buffer)==FALSE | sea_mask==TRUE){
    cat("Preparing data and loading DEM...")
    cat("\n")
    ## dem is GMTED2010
    dem=NULL
    data(dem,envir = environment())
    if(is.null(raster_extent)) v=extent(xy) else v=extent(raster_extent)
    dem=crop(dem,v)
    dem1 <- projectRaster(dem, crs=proj4)
    #plot(dem1)
  }
  ###
  
  if(grepl(longlat,proj4)==FALSE){
    dat2=data.frame(project(xy, proj4))} else dat2=data.frame(xy)
  
  colnames(dat2)=c("x","y")
  if(is.null(raster_extent)) {
    e <- extent(c(round(range(dat2$x)*2)/2,round(range(dat2$y)*2)/2))
  } else {
    toto=cbind(x=c(raster_extent[c(1,1,2,2)], raster_extent[1]+(raster_extent[2]-raster_extent[1])/2),
               y=c(raster_extent[c(3,4,3,4)], raster_extent[3]))
    e <- extent(project(toto, proj4))
  }
  
  
  # Cell sizes
  if(is.null(cell_size)){
    if(is.null(cell_sizex)) cell_sizex=200000
    if(is.null(cell_sizey)) cell_sizey=200000
  } else {
    cell_sizex=cell_size
    cell_sizey=cell_size
  }
  # Number of rows and columns
  nc=ceiling((e@xmax-e@xmin)/cell_sizex)
  nr=ceiling((e@ymax-e@ymin)/cell_sizey)
  
  r <- raster(e, ncol=nc, nrow=nr, resolution=c(cell_sizex,cell_sizey)) # Empty raster
  projection(r)<-proj4
  
  if(is.null(distance_buffer)) distance_buffer=300000
  
  ## Elevation stuff (median elevation in each predicted cell)  
  if(is.null(elevation_buffer)==FALSE | sea_mask==TRUE){
    temp=rasterToPoints(dem1)
    temp1=rasterize(temp[, 1:2], r, temp[,3], fun=mean)
    ## ???????? Hack for Sea Mask NEW 12 may 14
    dem2=dem1
    dem2[is.na(dem2)]=-9999
    dem2=rasterToPoints(dem2)
    temp2=rasterize(dem2[, 1:2], r, dem2[,3], fun=median)
    ## !!!!!!!!
    z=raster::intersect(temp1,r)
    # plot(temp1)
    elev1=rasterToPoints(temp1)[,3]
    dat1=as.data.frame(rasterToPoints(z))
  } else dat1=as.data.frame(rasterToPoints(r))
  
  
  ## Initiate a percentage counter
  if(verbose==TRUE){
    percent=seq(10,100,by=10)
    values=round(percent*length(dat1[,1])/100)
    cat("Spatio-temporal weighted interpolation...")
    cat("\n")
    cat("Percentage done: ")
  }
  
  #ptm <- proc.time()
  
  ## Main loop time--distance weighting 
  dat1[,3]=c()
  for(i in 1:length(dat1[,1])){
    if(grepl(longlat,proj4)){
      d=pointDistance(dat1[i,1:2], dat2,longlat=TRUE)} else {
        d=pointDistance(dat1[i,1:2], dat2,longlat=FALSE)}
    
    d=cbind(dat2,data[,3],d,triCube(d,distance_buffer))
    
    ## Elevation range search
    if(is.null(elevation_buffer)==FALSE ){
      elev2=extract(dem1,d[,1:2])
      elev1=extract(temp1,dat1[i,1:2])
      elev2[elev2<elev1-elevation_buffer]=NA
      elev2[elev2>elev1+elevation_buffer]=NA
      d=cbind(d,elev2)
      d[is.na(d[,6]),5]=0
    }
    d1=d[d[,5]!=0,]
    
    ## Time weight
    d1=cbind(d1,triCube(age-d1[,3],time_buffer))
    # head(d)
    if(is.null(elevation_buffer)==FALSE){
      colnames(d1)=c("x","y","age","dist","dweight","elev","tweight")
    } else colnames(d1)=c("x","y","age","dist","dweight","tweight")
    
    d1$weight=d1$dweight*d1$tweight
    d1$data=data[d[,5]!=0,4]
    
    ## Combine time and distance weights
    threshold=0.5
    if(sum(d1$weight,na.rm=TRUE)>threshold){
      dat1[i,3]=sum(d1$data*d1$weight,na.rm=TRUE)/sum(d1$weight,na.rm=TRUE)
    } else dat1[i,3]=NA
    
    ## Verbose
    if(i %in% values & verbose==TRUE)
      cat(percent[values==i]," ",sep="")
    #cat(i," ")
  }
  
  r1=rasterize(dat1[, 1:2], r, dat1[,3], fun=mean) # Fill raster with mean value 
  #   plot(r1)
  
  
  ## SEA MASK  
  if(sea_mask==TRUE){
    r2 <- temp2<(-1000) ## sea is now 1
    r2[r2==1]=NA
    ## Other mask (e.g. ice) see help for details
    if(is.null(other_mask)==FALSE){
      #       plot(r2) 
      r3=mask(r2,other_mask)
      #       plot(r3)
      r3[r3==0]=1
      r3[is.na(r3)]=0
      #       plot(r2-r3)
      r2=r2-r3;r2[r2!=0]=NA
      #       plot(r2)
    }
    # !!!!!!MASK
    r1=(r1-r2)
  }
  
  # plot(r1-r2)
  
  dat1=as.data.frame(rasterToPoints(r1))
  
  ## End of DATA part --------------------------------------------------------------------
  
  out=list(raster=r1,df=dat1,proj4=proj4,extent=e,points=dat2,res=c(cell_sizex,cell_sizey))
  class(out)="pfGridding"
  return(out)
  
}





#' Plot a "pfGridding" object.
#' 
#' Plot maps presenting gridded and transformed charcoal values obtained from
#' the \code{\link{pfGridding}} function.
#' 
#' @method plot pfGridding
#' @export
#' @param x An object returned by \code{\link{pfGridding}}.
#' @param continuous Logical, plot continuous (TRUE) or discrete (FALSE) colors
#' on the map.
#' @param col_class Numeric, if continuous is false define here color classes
#' (single values: col_class=5, or sequences col_class=seq(-15,15,5) are
#' accepted.)
#' @param col_lim Numeric, limits for plotting grid cells values, grid cells
#' with values beyond col_lim are not plotted.
#' @param xlim Numeric, map limits.
#' @param ylim Numeric, map limits.
#' @param empty_space Percentage, define empty space around the map.
#' @param cpal String, color palette to use see
#' \code{\link[RColorBrewer]{brewer.pal}}
#' @param anomalies Logical, adapt output for plotting anomalies or not (color
#' classes, etc..)
#' @param file Path/Filename.tiff, the function can output a GeoTiff file if
#' desired.
#' @param points Logical, plot charcoal sites on the map?
#' @param add An object of the class "SpatialPolygonsDataFrame" (sp) to be
#' ploted on the map.
#' @param add_color Color of the added SpatialPolygonsDataFrame.
#' @param plot_countries Logical, default FALSE (if TRUE plot countries borderlines and coastlines) 
#' @param \dots \dots{}
#' @return A ggplot2 "gg" object that could be further manipulated.
#' @author O. Blarquez
#' @seealso \code{\link{pfGridding}}
#' @examples
#' 
#' ID=pfSiteSel(id_region=="ENA0", l12==1, long>-85)
#' 
#' TR=pfTransform(ID,method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,4000))
#' 
#' p=pfGridding(TR,age=1000)
#' 
#' plot(p,empty_space=100)
#' 
#' # require(ggplot2)
#' # pp=plot(p,empty_space=100)
#' # pp+ggtitle("my title..")
#' 
plot.pfGridding=function(x,continuous=TRUE,
                         col_class=NULL,
                         col_lim=NULL,
                         xlim=NULL,ylim=NULL,empty_space=10,
                         cpal="YlGn",
                         anomalies=TRUE,
                         file=NULL,points=FALSE,add=NULL,add_color="white",plot_countries=FALSE,...){
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    install.packages("RColorBrewer")
  }
  
  y=long=lat=group=res1=res2=NULL ##  no visible binding for global variable 'y' ?
  
  x$df=as.data.frame(rasterToPoints(x$raster))
  x$df[,3][is.infinite(x$df[,3])] = NA #define Inf as NA 

  # Define classes for colors
  if(is.null(col_class)){
    if(anomalies==TRUE){
      b1=seq(floor(min(x$df$layer,na.rm=TRUE)),0,len=5)
      b2=seq(0,ceiling(max(x$df$layer,na.rm=TRUE)),len=5)
      breaks=c(b1,b2[2:5])
      breaks=unique(breaks)
    } else  breaks=seq(floor(min(x$df$layer,na.rm=TRUE)),ceiling(max(x$df$layer,na.rm=TRUE)),len=10)
  }
  if(is.numeric(col_class) & length(col_class)==1){
    if(anomalies==TRUE){
      b2=seq(0,ceiling(max(abs(x$df$layer)))+col_class,by=col_class)
      breaks=c(-rev(b2),b2[2:length(b2)])
    } else breaks=seq(floor(min(x$df$layer,na.rm=TRUE)),ceiling(max(x$df$layer,na.rm=TRUE))+col_class,by=col_class)
  }
  if(is.numeric(col_class) & length(col_class)>1){
    breaks=col_class
  }
  # End of options
  
  ## Define color limils
  if(is.null(col_lim))
    col_lim=c(floor(min(x$df$layer,na.rm=TRUE)),ceiling(max(x$df$layer,na.rm=TRUE)))
  
  x$df=cbind(x$df,class=cut(x$df$layer,breaks))
  #   data=na.omit(data)
  
  coast=NULL ## Load coast data and transform
  data(coast,envir = environment())
  coast=coast[coast$X<=180,]
  # plot((coast$X),(coast$Y),type="l")  
  xy=cbind(coast[,2],coast[,1])
  longlat="+proj=longlat|+proj = longlat|+proj =longlat|+proj= longlat|+proj=longlat"
  if(grepl(longlat,x$proj4)==FALSE){
    coast=data.frame(project(xy, x$proj4))} else coast=data.frame(xy)
  
  
  colnames(coast)=c("x","y")
  
  
  countries=NULL  ## Load countries data and transform
  data(countries,envir = environment())
  # plot((coast$X),(coast$Y),type="l")  
  xy=cbind(countries$x,countries$y)
  xy=xy[xy[,1]<=180,]
  countries=data.frame(project(xy, x$proj4))
  colnames(countries)=c("x","y")
  
  
  ## LIMITS
  if(is.null(xlim)){
    xplus=(x$extent@xmax-x$extent@xmin)*empty_space/100
    xlim=c(x$extent@xmin-xplus,x$extent@xmax+xplus)}
  
  if(is.null(ylim)){
    yplus=(x$extent@ymax-x$extent@ymin)*empty_space/100
    ylim=c(x$extent@ymin-yplus,x$extent@ymax+yplus)}
  
  ## Crop coast using limits
  coast=coast[coast$x>xlim[1]-8000000 & coast$x<xlim[2]+8000000 &
                coast$y>ylim[1]-8000000 & coast$y<ylim[2]+8000000,]
  #plot(coast[,1],coast[,2],type="l")
  
  x$points=data.frame(na.omit(x$points))
  colnames(x$points)=c("x","y")
  
  if(anomalies==TRUE) {
    if(cpal=="YlGn") {cpal="RdBu"}
    pal=rev(RColorBrewer::brewer.pal(9,cpal))} else pal=RColorBrewer::brewer.pal(9,cpal)
  
  testcol  <- colorRampPalette(pal)
  
  #   ## LIMITS
  #   if(is.null(xlim)){
  #     xplus=(x$extent@xmax-x$extent@xmin)*0.1
  #     xlim=c(x$extent@xmin-xplus,x$extent@xmax+xplus)}
  #   
  #   if(is.null(ylim)){
  #     yplus=(x$extent@ymax-x$extent@ymin)*0.1
  #     ylim=c(x$extent@ymin-yplus,x$extent@ymax+yplus)}
  
  ## SAME COLORS 
  if(continuous==FALSE & is.numeric(col_class) & length(col_class)>1){
    options(warn=-1)
    c1=breaks+(diff(breaks)/2)
    c3=cbind(xlim[2]+1e+6,ylim[2]+1e+6,c1)
    c3=c3[1:(length(c3[,1])-1),]
    c3=data.frame(c3,rr=(cut(c3[,3],breaks)))
    colnames(c3)=colnames(x$df)
    x$df=rbind(x$df,c3)
  }
  ## 
  pale=testcol(length(levels(x$df$class)))
  if(cpal=="Greys") {pale=rev(testcol(length(levels(x$df$class))+1)[1:(length(levels(x$df$class)))])}
  
  #display.brewer.pal(12,"Spectral")
  #pal=c(pal[9],pal[8],pal[6:1])
  ## On fait une carte avec ggplot2
  #pdf(file="/Users/Olivier/Desktop/x$dfMap.pdf",height=6,width=9)
  
  ## Add shp to the plot
  if(is.null(add)==FALSE){
    theclass=lapply(add@data, class) 
    theclass=which(theclass=="factor")
    add@data$id = add@data[,as.numeric(theclass[1])]
    add.points = fortify(add, region="id")
    add.df = plyr::join(add.points, add@data, by="id")    
  }
  
  ##
  x$df$res1=x$res[1]; x$df$res2=x$res[2]
  
  ## Plot only points
  if(points=="Only"){
    p=ggplot(x$df) +
      geom_polygon(data=coast,aes(x=x,y=y),colour="grey80",fill="grey80")+
      coord_cartesian(xlim=xlim,ylim=ylim)+xlab("Longitude")+ylab("Latitude")+
      theme_bw(base_size = 16)+
      geom_point(data=x$points,aes(x=x,y=y),colour="grey40")
    if(is.null(add)==FALSE){
      p=p+geom_polygon(data=add.df,aes(x=long,y=lat,group=group),fill=add_color,colour="grey80")
    }
    if (plot_countries==TRUE) {
      p=p+geom_path(data=countries,aes(x=x,y=y),colour="white")
    }
  } else {
    ## Plot interp data
    if(continuous==FALSE){
      p=ggplot(x$df) 
      if(cpal=="Greys") {p=p+geom_polygon(data=coast,aes(x=x,y=y),colour="black",fill="white")
      } else {p=p+geom_polygon(data=coast,aes(x=x,y=y),colour="grey80",fill="grey80")}
      p=p+geom_tile(data=x$df,aes(x=x, y=y, fill = class,width=res1,height=res2))+
        scale_fill_manual(values = pale,name="")+
        coord_cartesian(xlim=xlim,ylim=ylim)+xlab("Longitude")+ylab("Latitude")+
        theme_bw(base_size = 16)
      if(points==TRUE) p=p+geom_point(data=x$points,aes(x=x,y=y),colour="grey40")
      if(is.null(add)==FALSE){
        p=p+geom_polygon(data=add.df,aes(x=long,y=lat,group=group),fill=add_color,colour="grey80")
      }
      if (plot_countries==TRUE) {
        p=p+geom_path(data=countries,aes(x=x,y=y),colour="white")
      }
    } else {
      p=ggplot(x$df) +
        geom_polygon(data=coast,aes(x=x,y=y),colour="grey80",fill="grey80")+
        geom_tile(data=x$df,aes(x=x, y=y, fill = layer, width=res1, height=res2))+
        scale_fill_gradient2(high=pal[9],low=pal[1],mid="white",limits=col_lim)+
        coord_cartesian(xlim=xlim,ylim=ylim)+xlab("Longitude")+ylab("Latitude")+
        theme_bw(base_size = 16)
      if(points==TRUE) p=p+geom_point(data=x$points,aes(x=x,y=y),colour="grey40")
      if(is.null(add)==FALSE){
        p=p+geom_polygon(data=add.df,aes(x=long,y=lat,group=group),fill=add_color,colour="grey80")
      }
      if (plot_countries==TRUE) {
        p=p+geom_path(data=countries,aes(x=x,y=y),colour="white")
      }
    }
  }
  
  p
  if(is.null(file)==FALSE){
    projection(x$raster)<-sp::CRS(x$proj4)
    writeRaster(x$raster, filename=file, format="GTiff",overwrite=TRUE)
  }
  return(p)
}

## -------------------------------------------------------------------------------------------




#' Tukey's Tricube weight function
#' 
#' From the EGRET package http://usgs-r.github.io/EGRET/ Robert Hirsch and
#' Laura De Cicco
#' 
#' Computes the tricube weight function on a vector of distance values (d),
#' based on a half-window width of h, and returns a vector of weights that
#' range from zero to 1.
#' 
#' @param d numeric vector of distances from the point of estimation to the
#' given sample value
#' @param h numeric value, the half-window width, measured in the same units as
#' d
#' @return w numeric vector of weights, all 0<=w<=1
#' @keywords statistics weighting
#' @examples
#' 
#' h<-10
#'  d<-c(-11,-10,-5,-1,-0.01,0,5,9.9,10,20)
#'  triCube(d,h)
#' 
triCube<-function(d,h) {
  #  triCube is Tukey tricubed weight function
  #    first argument, d, is a vector of the distances between the observations and the estimation point
  #    second argument, h, is the half window width
  #    it returns a vector of weights (w) for the observations in the vector, d
  n<-length(d)
  zero<-rep(0,n)
  ad<-abs(d)
  w<-(1-(ad/h)^3)^3
  w<-pmax(zero,w)
  return(w)
}
