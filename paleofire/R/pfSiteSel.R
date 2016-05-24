#' GCD sites selection methods
#' 
#' Main function used for site selection, uses data strored in
#' data(paleofiresites) to perform site selection according to multiple
#' criterion, those criterions could be either geographic, based on series
#' attributes (e.g. # of datings), or on sites attributes (e.g. biome).
#' 
#' Use data(paleofiresites);names(paleofiresites) to retrieve the conditions
#' that could be used to select sites i.e.: id_site, site_name, lat, long,
#' elev, pref_units, biome, id_region, id_country, id_site_type, water_depth,
#' id_basin_size, id_catch_size, id_land_desc, dating_type, min_est_age,
#' max_est_age, num_dating, age_model, data_source, qtype, rf99, l12, num_samp,
#' date_int.
#' 
#' @param ... Any combination of conditions defined by relational operators and
#' or logical operators that are applied on the "paleofiresites" dataset. See
#' examples below:
#' @return An object of the class "pfSiteSel" (list) with "id_site" and
#' "site_name" components.
#' @author O. Blarquez
#' @seealso \code{\link[GCD]{paleofiresites}}
#' @examples
#' 
#' ## Sites selection examples
#' 
#' ## Select all sites
#' ID=pfSiteSel()
#' 
#' ## Site in the Biome #8
#' ID=pfSiteSel(biome==8)
#' plot(ID,zoom="world")
#' 
#' ## Site in the biome #8 or in the biome #6
#' ID=pfSiteSel(biome==8 | biome==6)
#' 
#' ## Sites in North America by geographic location
#' ID=pfSiteSel(lat>25, lat<75, long<(-45), long>-150) 
#' plot(ID,zoom="world")
#' 
#' ## is equivalent to:
#' ID=pfSiteSel(lat>25 & lat<75 & long<(-45) & long>-150) 
#' plot(ID,zoom="world")
#' 
#' ## By region criterion
#' ID=pfSiteSel(id_region==c("ENA0","WNA0"))
#' plot(ID,zoom="world")
#' 
#' ## WRONG, use the %in% operator when concatenating two characters
#' # ID=pfSiteSel(id_region %in% c("ENA0","WNA0"))
#' # plot(ID,zoom="world")
#' 
#' ## Pas-de-Fond site
#' pfSiteSel(site_name=="Pas-de-Fond")
#' 
#' ## All sites in  eastern North America that are not Pas-de-Fond
#' pfSiteSel(site_name!="Pas-de-Fond", id_region=="ENA0")
#' 
#' ## Sites with on average one dating point every 250 to 300 yrs
#' pfSiteSel(date_int>=250 & date_int<=300)
#' 
#' ## Sites between 0, 100 m elevation in Asia
#' ID=pfSiteSel(elev>0 & elev<100, id_region=="ASIA")
#' 
#' ## All sites that are not marine nor fluvial
#' ID=pfSiteSel(id_land_desc!="MARI" , id_site_type!="FLUV" & id_site_type!="LFLU")
#' plot(ID)
#' 
#' 
pfSiteSel <- function(...) {
  
  ## Load data (bindind...)
  paleofiresites<-NULL
  data(paleofiresites,envir = environment())
  
  ## The eval function:
  theeval=function(thelist){
    eval(thelist,paleofiresites)
  }
  
  ## Retrieve all args in dots
  args=eval(substitute(alist(...)))
  
  if(length(args)>0){
    ## Eval all args
    c=lapply(args,theeval)
    d=matrix(unlist(c),ncol=length(args))
    d[is.na(d)]=FALSE
    
    ## All TRUE?
    finalTF=ifelse(rowSums(d == TRUE) == length(args), TRUE, FALSE)
    id=paleofiresites[finalTF,]$id_site
  } else id=paleofiresites$id_site
  
  ## Output:
  site_name=as.character(paleofiresites$site_name[paleofiresites$id_site %in% id])
  
  output=list(id_site=id,site_name=site_name)
  class(output)="pfSiteSel"
  return(output)
}

## Summary function




#' summary.pfSiteSel
#' 
#' Return a summary table for an object of the class "pfSiteSel"
#' 
#' @method summary pfSiteSel
#' @export
#' @param object An object of the class "pfSiteSel".
#' @param \dots \dots{}
#' @return Data.frame, returns the following informations: "id_site", "lat",
#' "long" "elev", "min_est_age", "max_est_age", "num_dating", "date_int",
#' "num_samp", "l12", "rf99".
#' @author O. Blarquez
#' @examples
#' 
#' ID=pfSiteSel(id_site==2)
#' summary(ID)
#' 
summary.pfSiteSel=function(object,...){
  
  ## Avoid no visible binding for global variable
  paleofiresites=NULL; rm(paleofiresites)
  paleofiredata=NULL; rm(paleofiredata)
  
  data(paleofiresites,envir = environment())
  data(paleofiredata,envir = environment())
  coln=length(paleofiresites[1,])
  
  table=paleofiresites[paleofiresites$id_site %in% object$id_site,]
  rownames(table)=table$site_name
  table=subset(table, select=c("id_site", "lat", "long",
                               "elev", "min_est_age", "max_est_age", 
                               "num_dating",  "date_int", "num_samp", "l12", "rf99"))
  #print(table)
  return(table)
}

## Plot functions




#' plot.pfSiteSel
#' 
#' Plot an object of the class "pfSiteSel"
#' 
#' @method plot pfSiteSel
#' @export
#' @param x An object of the class "pfSiteSel".
#' @param add An object returned by pfAddData (optional).
#' @param type Character, type of plot among "Map" or "Chronology".
#' @param zoom Character, zooming factor for type="Map": "Sites" or "World"
#' @param pch Pointer type see \code{\link[graphics]{plot}}.
#' @param xlim Numeric, x axis limits.
#' @param ylim Numeric, y axis limits.
#' @param cex Numeric, size of points.
#' @param plot_countries Logical, default FALSE (if TRUE plot countries borderlines and coastlines) 
#' @param \dots \dots{}
#' @author O. Blarquez
#' @examples
#' 
#' ID=pfSiteSel(id_region=="ENA0", long>-100)
#' plot(ID,zoom="world")
#' 
#' 
plot.pfSiteSel=function(x,add=NULL,type="Map",zoom="Sites",pch="|",
                        xlim=NULL, ylim=NULL, cex=1, plot_countries=FALSE,...)
  
{
  ## Avoid no visible binding for global variable
  paleofiresites=NULL; rm(paleofiresites)
  coast=countries=NULL; rm(coast,countries)
  
  data(paleofiresites,envir = environment())
  data(coast,envir = environment())
  data(countries,envir = environment())
  
  ## Chronology
  if(type=="Chronology"){
    data(paleofiredata,envir = environment())
    paleofiredata=paleofiredata[paleofiredata$id_site %in% x$id_site,]
    
    IDsorted=data.frame(IDs = c(x$id_site),
                        Lat = c(paleofiresites[paleofiresites$id_site %in% x$id_site,]$lat),
                        labels = as.character(paleofiresites$site_name[paleofiresites$id_site %in% x$id_site]))
    
    IDsorted=IDsorted[with(IDsorted, order(Lat)), ]
    ## Xlim
    if(is.null(xlim)) xlim=range(paleofiredata$EST_AGE)
    
    ## Plot
    par(mar=c(4,14,2,8))
    plot(NULL, type = "n", 
         ylim = c(1,length(x$id_site)),xlim=xlim,axes=FALSE,ylab="",xlab="Age",main="Sampling resolution")
    n=c()
    for(i in 1:length(x$id_site)){
      samples=paleofiredata$EST_AGE[paleofiredata$id_site %in% IDsorted$IDs[i]]
      points(samples,rep(i,length(samples)),pch=pch,cex=cex)
      n[i]=length(samples)
    }
    axis(2, at=seq(1,length(IDsorted$IDs),1), labels = FALSE)   
    IDsorted$labels=gsub("[\x87]", "c", IDsorted$labels) 
    IDsorted$labels=gsub("[\x85]", "a", IDsorted$labels) 
    IDsorted$labels=gsub("[\x82]", "e", IDsorted$labels) 
    IDsorted$labels=gsub("[\x8a]", "e", IDsorted$labels) 
    text(y = seq(1,length(IDsorted$IDs),1), par("usr")[1], labels = IDsorted$labels, srt = 0, pos = 2, xpd = TRUE)
    axis(side = 1, at = seq(0, 99000, by = 500), labels = FALSE, tcl = -0.2) 
    axis(4, at=seq(1,length(IDsorted$IDs),1), labels = FALSE)    
    text(y = seq(1,length(IDsorted$IDs),1), par("usr")[2], labels = paste(round(IDsorted$Lat,digits=1),"/",n,sep=""), srt = 0, pos = 4, xpd = TRUE)
    
    paste(round(IDsorted$Lat,digits=1),"/",n,sep="")
    axis(1)
  }
  
  ## MAPS
  if(type=="Map"){
    
    if(zoom=="World"|zoom=="world"){
      plot(paleofiresites$long,paleofiresites$lat,
           col="blue",xlab="Longitude",ylab="Latitude")
      
     if (plot_countries==TRUE) {
       lines(countries$x,countries$y)
     } else lines(coast$X,coast$Y)
      
      points(paleofiresites[paleofiresites$id_site %in% x$id_site,]$long,
             paleofiresites[paleofiresites$id_site %in% x$id_site,]$lat, 
             bg="red",col = "red",pch = 21,xlab="Longitude",ylab="Latitude")
      if(is.null(add)==FALSE)
        points(add$metadata$LONGITUDE,
               add$metadata$LATITUDE, 
               bg="red",col = "red",pch = 21)
    }
    
    if(zoom=="Sites"|zoom=="sites"){
      # Draw map
      xl=as.vector(paleofiresites[paleofiresites$id_site %in% x$id_site,]$long)
      yl=as.vector(paleofiresites[paleofiresites$id_site %in% x$id_site,]$lat)
      
      if(is.null(xlim))
        xlim=range(xl[!is.na(xl) & is.finite(xl)])
      if(is.null(ylim))
        ylim=range(yl[!is.na(yl) & is.finite(yl)])
      
      
      plot(paleofiresites$long,paleofiresites$lat,col="blue",xlab="Longitude",ylab="Latitude",xlim=xlim,ylim=ylim)
      points(paleofiresites[paleofiresites$id_site %in% x$id_site,]$long,paleofiresites[paleofiresites$id_site %in% x$id_site,]$lat,bg="red",col = "red",pch = 21)
      
      if (plot_countries==TRUE) {
        lines(countries$x,countries$y)
      } else lines(coast$X,coast$Y)
      
      if(is.null(add)==FALSE)
        points(add$metadata$LONGITUDE,
               add$metadata$LATITUDE, 
               bg="red",col = "red",pch = 21)
      
    }
  }
}

