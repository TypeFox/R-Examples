# Purpose        : Convert SGDF/raster map to a polygon map (grid cells);
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : Dylan Beaudette (debeaudette@ucdavis.edu); Pierre Roudier (pierre.roudier@landcare.nz); 
# Status         : working version
# Note           : Not recommended for large grids;

grid2poly <- function(obj, var.name = names(obj)[1], reproject = TRUE, method = c("sp", "raster", "RSAGA")[1], tmp.file = TRUE, saga_lib = "shapes_grid", saga_module = 3, silent = FALSE, ...){

    # print warning:
    if(length(obj)>1e4){
    warning("Operation not recommended for large grids (>>1e4 pixels).", immediate. = TRUE)
    }
       
    if(method=="raster"){
        r <- raster(obj[var.name])
        pol <- rasterToPolygons(r)
        names(pol) <- var.name
    }
    
    else{
    if(method=="RSAGA"){
      if(!rsaga.env()[["cmd"]]=="NULL"){
        if(tmp.file==TRUE){
          tf <- tempfile() 
        } else { 
          tf <- var.name
        }

        # first, write SGDF to a file:
        obj <- as(obj[var.name], "SpatialPixelsDataFrame")
        writeGDAL(obj[var.name], paste(tf, ".sdat", sep=""), "SAGA")
        # saga_lib name and saga_module might change in the future versions of SAGA!
        # SAGA GIS 2.0.8
        rsaga.geoprocessor(lib=saga_lib, module=saga_module, param=list(GRIDS=paste(tf, ".sgrd", sep=""), SHAPES=paste(tf, ".shp", sep=""), NODATA=TRUE, TYPE=1), show.output.on.console = silent)
        if(requireNamespace("maptools", quietly = TRUE)){
          pol <- maptools::readShapePoly(paste(tf, ".shp", sep=""), proj4string=obj@proj4string)
        } else {
          pol <- readOGR(paste(tf, ".shp", sep=""))
        }
      }
        
        else { stop("SAGA GIS path could not be located. See 'rsaga.env()' for more info.") }
    }
    
        else {
        obj <- as(obj[var.name], "SpatialPixelsDataFrame")
        # pol <- as.SpatialPolygons.SpatialPixels(obj) # EJP: deprecated
        # pol <- SpatialPolygonsDataFrame(pol, data=data.frame(var.name = obj@data[,var.name]), match.ID=FALSE)
		    pol = as(obj, "SpatialPolygonsDataFrame")
    }
    }
    
    # Checking projection:
    prj.check <- check_projection(pol, control = TRUE)

    # Trying to reproject data if the check was not successful
    if (!prj.check&reproject==TRUE) {  pol <- reproject(pol)  }

    # Convert to SPolyDF:
    dm <- data.frame(obj@data[,var.name])
    names(dm) <- var.name
    pol <- SpatialPolygonsDataFrame(pol, dm, match.ID=FALSE)
    
    return(pol)
} 

# end of script;
