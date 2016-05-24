RasterizeTraj <-
function(spLines, resolution=10000, reduce = TRUE, parallel=FALSE) {
    # This function produces a grid over an specified area and then computes the 
    # frequency of lines that cross the grid cells. 
    # 
    # Args:
    #   spLines: An object of class Spatial Lines created by the function Df2SpLines
    #   reduce: Boolean: If TRUE the result will be reduced to one raster object
    #           if FALSE, this function will return a list of raster Layers. The size of 
    #           the list is equal to the number of available cores in the system
    #   
    # Returns:
    #   An obejct of class RasterLayer
        
    
    getRasterGrid <- function(sp.lines, xmn, xmx, 
                              ymn, ymx, ncols=40,
                              nrows=40, resolution=10000, ext=ext)
    {
      # create raster object
      rast <- raster(xmn = xmn, 
                     xmx = xmx, 
                     ymn = ymn, 
                     ymx = ymx, 
                     ncols = ncols, 
                     nrows = nrows)
      
      
      # set all the grids to NA
      rast <- setValues(rast, NA)
      
      crs1 <- "+proj=aea +lat_1=46 +lat_2=60 +lat_0=44 +lon_0=-68.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
      rast <- projectRaster(rast, crs = crs1, res = resolution)
      
      # get the projection from the sp object
      crs2 <- proj4string(sp.lines)
      
      # Reproject
      rast <- projectRaster(rast, crs = crs2)
      
      rast
    }
    
    # get the bounding box of the spLines object
    ext <- extent(spLines)
    
    if(parallel == TRUE) {
    
      # split the sp lines object into N set of lines
      # where N is the number of cores available
      cores <- detectCores()
      
      list.splines <- SplitSpLines( spLines, cores )
      
      cl <- makeCluster(cores)
      
      registerDoParallel(cl)
      
      rast2 <- foreach(sp.lines = list.splines, .combine='c', .packages="raster") %dopar%
      {
        # And then ... rasterize it! This creates a grid version 
        # of your points using the cells of rast, values from the IP field:
        rast <- getRasterGrid(sp.lines, 
                              xmn=xmin(ext),
                              xmx = xmax(ext),
                              ymn = ymin(ext),
                              ymx = ymax(ext),
                              resolution=resolution)
        
        rasterize(sp.lines, rast, fun='count', background=0) 
      }
      
      stopCluster(cl)
      
      if(reduce == T){
        rast2 <- Reduce("+", rast2)
        
        # replace all 0 values per NA
        rast2[rast2==0] <- NA
        
      }
      
      rast2
      
    } else {
        # And then ... rasterize it! This creates a grid version 
        # of your points using the cells of rast, values from the IP field:
        rast <- getRasterGrid(spLines, 
                              xmn=xmin(ext),
                              xmx = xmax(ext),
                              ymn = ymin(ext),
                              ymx = ymax(ext),
                              resolution=resolution)
        
        rast2 <- rasterize(spLines, rast, fun='count') 
        
        rast2
    }

}
