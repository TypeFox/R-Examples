# Remarks
# Is it reasonable to have ainfo <<- read.area.info(finfo,...) to make sure that ainfo 
# is also available at the top level after being delivered to function read.areas?

readAreaInfo = function(fname = "ainfo.txt", id = "id", 
            iobs = "iobs",obs="obs",unc="unc",
            filenames="filenames", sep = "\t", debug.level = 1, 
            moreCols = list(NULL)) {
# Separate function to read in information about the areas, with possibility to define column names
# fname = name of file with information
# inum = area number
# id = internal identity - the way areas are stored on hard drive
# iobs = number of observations for each station
# obs = the actual observation
# unc = The uncertainty of the observation, standard deviation
#       This variable is optional
# filenames = filenames for areas
# MoreCols = other variables the user wants to pass on to ainfo
#  cat(paste(fname))
  if (debug.level > 1)  print(paste(fname, id, iobs, obs, unc, 
                       filenames, debug.level, moreCols))
  ainfot <- read.csv(fname,header=TRUE,sep = sep)
  if (debug.level > 1) print(summary(ainfot))
  ainfo = data.frame(id = ainfot[,names(ainfot) == id],
                     iobs = ainfot[,names(ainfot) == iobs],
                     obs = ainfot[,names(ainfot) == obs])
  if (unc %in% names(ainfot) ) 
      ainfo = data.frame(ainfo, unc = ainfot[,names(ainfot) == unc])
  if (filenames %in% names(ainfot) ) 
      ainfo = data.frame(ainfo, 
            filenames = ainfot[,names(ainfot) == filenames])
# Including remaining arguments, if user wants to include
  ncols = length(moreCols)
  if (!is.null(moreCols[[1]])) {
    for (icol in 1:ncols) {
      col = moreCols[[icol]]
      if(debug.level >1) cat(paste("moreCols",icol,col,"\n"))
      ainfo = data.frame(ainfo,col1 = ainfot[,names(ainfot)==col])
      names(ainfo) = c(names(ainfo[1:(length(names(ainfo))-1)]),col)
      if(debug.level >1) cat(paste("moreCols2",icol,col,"\n"))
      if(debug.level >1) cat(paste("names(ainfo)",names(ainfo),"\n"))
    }
  }
  return(ainfo)
}




readAreas = function(object, adir=".",ftype = "xy",projection = NA, ...)  {
# ainfo is 1 - ainfo e.g. read by readAreaInfo
#          2 - name of the file to pass to readAreaInfo. ainfo is in that case delivered as a top level data.frame
# pdif gives directory to areal information
# need option to use other separators as well
# Output of this function is a list consisting of
#      1 Updated version of ainfo
#      2 A list of Spatial polygons defining the borders of the areas or a set of Spatial grids
  if (is.character(object)) {
    cat(paste("calling readAreaInfo with filename ",object,"\n"))
    ainfo = readAreaInfo(object,...)
  } else {
    ainfo = object
  }
  cat(paste(names(ainfo)))
  cat(paste("\n"))
  if (sum(names(ainfo) == "filenames") == 1) {
    fnames = paste(adir,"/",ainfo$filenames,sep="")
  } else {
    fnames = ainfo$id
  }
  areas = list()
  row.names(ainfo) = c(1:dim(ainfo)[1])
  if (ftype == "xy") {
    fnames = paste(adir,"/",fnames,".xy",sep="")
    Srl = list()
    for (i in 1:length(fnames)) {
      cat(paste("reading first polygon",i,length(fnames),"\n"))
      boun = read.table(fnames[i],header = FALSE)
      names(boun) = c("x","y")
      coordinates(boun) = ~x+y
      boun = Polygon(boun)
      cat(paste("adding data to ainfo",i,boun@area,boun@labpt[1],"\n"))
      cat(paste(" Finished polygon\n"))
      Srl[[i]] = Polygons(list(boun),ID = as.character(i))
    }
    Sr = SpatialPolygons(Srl, proj4string=CRS(as.character(projection)))
#  } else if (ftype == "shp") {
      # This part is when each file has a single shape
      # NOT PROPERLY IMPLEMENTED - need testing with real shapes
      # The files will probably be read as lists of polygons, necessary
      # to extract the actual polygon
#    require(maptools)
#    fnames = paste(adir,"/",fnames,".shp",sep="")
#    Srl = list()
#    for (i in 1:dim(fnames)) {
#      boun = readShapePoly(fnames[i])
#      Srl[[i]] = Polygons(boun,ID = as.character(i))
#    }
#    Sr = SpatialPolygons(Srl, proj4string=CRS(as.character(projection)))
#  } else if (ftype == "shps") {
      #This clause is when one shapefile includes all the shapes
      # NOT properly tested yet
      # Necessary to split all shapes into single Polygons
#    require(maptools)
#    Sr = readShapePoly(adir)
#    SPDF = SpatialPolygonsDataFrame(Sr,data = ainfo)
  } else stop(paste("Filetype", ftype, "not recognized"))
  ainfo$area = unlist(lapply(Sr@polygons,FUN = function(poly) poly@area))
  ainfo$labx = unlist(lapply(Sr@polygons,FUN = function(poly) poly@labpt[1]))
  ainfo$laby = unlist(lapply(Sr@polygons,FUN = function(poly) poly@labpt[2]))
  ainfo$bdim = unlist(lapply(Sr@polygons,FUN = function(poly) dim(poly)[1]))
  SPDF = SpatialPolygonsDataFrame(Sr, data = ainfo, match.ID = TRUE)
  SPDF
}



