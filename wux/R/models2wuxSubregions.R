
# ----------------------------------------------------------------
# $Author: geh $
# $Date: 2015-03-31 14:00:14 +0200 (Tue, 31 Mar 2015) $
# $Rev: 331 $
# ----------------------------------------------------------------

## ----------------------------------------------------------------
## WUX subregion functions will do following homeworks:
##
## - Basically getting offset and count for NetCDF files.
## - Plus getting a mask array, indicating which values are NA
##     (i.e. ouside subregion).
## - Helperfunctions for rectangular clipping, shapefiles and
##      subregionfiles.
##
## Remark: Rectangles and shapefiles will
##         be coerced to a Polygon in the first place. After that
##         we check which data- or gridpoints lie inside this Polygon.
##         Subregionsfiles will be interpolated to the model grid
##         via nearest neighbor.
## ----------------------------------------------------------------

GetSubregionShapes <- function(subregions,
                               grid.filenames,
                               area.fraction,
                               spatial.weighting,
                               lonlat.var.name, ...) {
  ## Get offset and count for reading in rectangular array from NetCDF
  ## file and index mask for clipping subregion from this rectangule.
  ## Rectangular subregions and shapefiles will be coerced to a polygon first,
  ## (call GetSubregionPolygons) and subregionsfiles will be interpolate
  ## directly to the model grid using nearest neighbor
  ## (call SubregionFile2ModelGrid).
  ##
  ## Args:
  ##   subregions: subregions tag from user.input.
  ##   grid.filenames: Filename of grid for specific climate model for which
  ##                   to get offset and count.
  ##
  ## Returns:
  ##   List of subregions containing lists with 'offset', 'count' and 'mask'
  ##   for exact clipping of subregion.
  ##
  ## History:
  ##   2010-11-24 | Original code (thm)
  ##   2011-02-14 | wrapping function (180 and 360 degrees) added (thm)
  ##   2011-06-01 | small bugfix for unnamed subregions in shapefiles (thm)
  ##   2011-07-08 | masking cells with area fraction of overlaying polygon
  ##                (mir & thm)
  ##
  ## TODO:
  ##   clean implementation of wrap.to180 objects...

  ## read in constant fields for climate model grid
  model.grid <- ReadGridFile(grid.filenames,
                             lonlat.var.name = lonlat.var.name,
                             ...)

### check what type of subregion. if rectangular then...,
  ## if shape then... if false then ..GetIndicesFromSubregiFile as usual

  ## initialize list containing clipping informations being returned
  output.list <- list()

  ## count over all subregions...
  for (subregion.counter in seq(along = subregions)){
    ## get subregion list entry
    subregion <- subregions[[subregion.counter]]
    ## get subregion name
    subregion.name <- names(subregions[subregion.counter])
    ## get type of subregion
    subregion.type <- NULL

    ## should be extra function

    ## Western, Eastern, Northern and Southern extensions
    if (length(subregion) == 4 && is.numeric(subregion))
      subregion.type <- "rectangle"
    ## shapefile
    if (length(subregion) %in% c(2,3,4) && is.list(subregion))
      subregion.type <- "shapefile"
    ## single subregionfile for all models. interpolating will be necessary...
    ## right now, we dont use that.
    if (length(subregion) >= 5 && is.list(subregion))
      subregion.type <- "single.subregionfile"
    ## subregionfile for every model in its subregiontag....
    if (is.null(subregion))
      subregion.type <- "multiple.subregionfiles"
    ##
    if (class(subregion) == "SpatialPolygons")
      subregion.type <- "SpatialPolygons"
    ## subregions tag in general.input has unknown form
    if (is.null(subregion.type))
      stop("UNKNOWN SUBREGION TYPE, CHECK IN GetSubregionPolygons")

### Make polygon, then get offset, count and clipping mask
###  for the polygon on the model grid.
    if (subregion.type %in% c("shapefile", "rectangle")){
      cat("\n      getting polygons for subregion \"", subregion.name,
          "\" (", subregion.type, ")\n", sep="")

      ## get polygons
      poly.list <- GetSubregionPolygons(subregion, subregion.type, ...)
      ## look what subregions longitudes to wrap to 360 or 180
      ## wrap all to 180 by default
      if ("subregionnames" %in% names(subregion)){
        wrap.to <- rep("180", length(subregion$subregionnames))
        names(wrap.to) <- subregion$subregionnames
      } else {
        wrap.to <- rep("180", length(poly.list))
      }
      if ("wrap.to" %in% names(subregion))
        wrap.to[names(subregion$wrap.to)] <- subregion$wrap.to

      ## getting constant fields of model
      cat("      get offset, count and mask for polygon(s) of \"",
          subregion.name, "\" on model grid\n", sep="")
      ## get indices for modelgrid from polygons
      clippings.subreg.list <- vector("list", length(poly.list)) # initialize
      names(clippings.subreg.list) <- names(poly.list)

      for (ii in seq(along=poly.list)) {
        sp.subregion <- poly.list[[ii]]
        wrap.to.subregion <- wrap.to[ii]

        ## wrap lon values to 0 to 360 degrees or -180 to 180 degrees
        if (any(!as.character(wrap.to) %in% c("180", "360")))
          stop("wrap.to NOT KNOWN. YOU HAVE TO TYPE 180 or 360.")
        if (wrap.to.subregion == "180"){
          sp.subregion <- wrapTo180(sp.subregion)
          model.grid$lon <- wrapTo180(model.grid$lon)
        }
        if (wrap.to.subregion == "360"){
          sp.subregion <- wrapTo360(sp.subregion)
          model.grid$lon <- wrapTo360(model.grid$lon)
        }

        ## if area.fraction is true, calculate the weighting mask
        ## else just check which cell centroids are inside shape
        if ( area.fraction ) {
          clippings.subreg.list[[ii]] <-
            GetWeightingMaskFromPolygon(sp.subregion, model.grid)
        } else {
          clippings.subreg.list[[ii]] <-
            GetIndicesFromPolygon(sp.subregion, model.grid)
        }

        ## if spatially weighted average is enabled, a weight matrix will be
        ## calculated using cos(lat). This is a valid approximation for
        ## regular grids only!
        if ( spatial.weighting ) {
          ## check if regular grid
          ## (longitudes have same distance for all latitudes)
          lon <- model.grid$lon
          lon.dist <- lon

          lat.dim <- seq(dim(lon)[2])   #2nd dim is latitude dim
          lon.dim <- seq(dim(lon)[1])   #1st dim is longitude dim
          for (i in seq(lon.dim))
            lon.dist[i, ] <- lon[max(c(i - 1, 1)), ] - lon[i, ]
          b <- NULL
          ## for each longitude (rows): all values equal?
          for (i in seq(lon.dim))
            b <- c(b, all(lon.dist[i, ] == (lon.dist[i, ][1])))
          is.regular.grid <- all(b)

          if ( is.regular.grid ) {
            ## get weighting cosinus weighting
            lat.radiant <- model.grid$lat * pi/180
            lat.weight <- cos(lat.radiant)
            ## cut out weightmatrix for corresponding domain
            offset <- clippings.subreg.list[[ii]]$offset.lon.lat
            count<- clippings.subreg.list[[ii]]$count.lon.lat
            mask <- clippings.subreg.list[[ii]]$mask
            lat.weight <- lat.weight[offset[1]:(offset[1] + count[1] - 1),
                                     offset[2]:(offset[2] + count[2]  - 1)]
            is.na(lat.weight) <- mask

            ## If there are no weights (e.g. area fraction weights),
            ## then use spatial weights as calculated here.
            ## Else multiply these weights.
            if ( is.null(clippings.subreg.list[[ii]]$weight) ) {
              clippings.subreg.list[[ii]]$weight <- lat.weight
            } else {
              clippings.subreg.list[[ii]]$weight <-
                clippings.subreg.list[[ii]]$weight * lat.weight
            }

          } else {
            stop( "ERROR: USING COSINE AREA WEIGHTING ON NON-REGULAR GRID!" )
          }
        }

        ## add lon lat info to list
        clippings.subreg.list[[ii]]$lon <-
          model.grid$lon[clippings.subreg.list[[ii]]$offset.lon.lat[1]:
                         (clippings.subreg.list[[ii]]$offset.lon.lat[1] +
                          clippings.subreg.list[[ii]]$count.lon.lat[1] - 1),
                         clippings.subreg.list[[ii]]$offset.lon.lat[2]:
                         (clippings.subreg.list[[ii]]$offset.lon.lat[2] +
                          clippings.subreg.list[[ii]]$count.lon.lat[2] - 1)]

        clippings.subreg.list[[ii]]$lat <-
          model.grid$lat[clippings.subreg.list[[ii]]$offset.lon.lat[1]:
                         (clippings.subreg.list[[ii]]$offset.lon.lat[1] +
                          clippings.subreg.list[[ii]]$count.lon.lat[1] - 1),
                         clippings.subreg.list[[ii]]$offset.lon.lat[2]:
                         (clippings.subreg.list[[ii]]$offset.lon.lat[2] +
                          clippings.subreg.list[[ii]]$count.lon.lat[2] - 1)]
      }

      ## NAME SUBREGIONS
      ##  In case subregions already have a name, concatenate it with subregion
      ##  name, else just give single subregion name.
      has.names <- !is.null(names(clippings.subreg.list))
      if (has.names){
        names(clippings.subreg.list) <- paste(subregion.name,
                                              names(clippings.subreg.list),
                                              sep=".")
      } else {
        names(clippings.subreg.list) <- subregion.name
      }
    }
### Interpolate with nearest neighbor to model grid, then
###  get offset, count and clipping mask.
    if (subregion.type %in% c("single.subregionfile")){
      cat("\n  interpolate subregions \"",
          subregion.name, "\" to model grid\n", sep="")
      ## get interpolated subregions to model grid via nearest neighbor
      subregions.list <- SubregionFile2ModelGrid(subregion, model.grid)
      ## get clipping indices for subregions
      cat("      get offset, count and mask for subregions \"",
          subregion.name, "\" on model grid\n", sep="")
      clippings.subreg.list <- lapply(subregions.list, GetOffsetCountAndMask,
                                      dims = model.grid$dims )

      ## add lon lat info to list
      for ( ii in 1:length(clippings.subreg.list) ) {
        clippings.subreg.list[[ii]]$lon <-
          model.grid$lon[clippings.subreg.list[[ii]]$offset.lon.lat[1]:
                         (clippings.subreg.list[[ii]]$offset.lon.lat[1] +
                          clippings.subreg.list[[ii]]$count.lon.lat[1] - 1),
                         clippings.subreg.list[[ii]]$offset.lon.lat[2]:
                         (clippings.subreg.list[[ii]]$offset.lon.lat[2] +
                          clippings.subreg.list[[ii]]$count.lon.lat[2] - 1)]

        clippings.subreg.list[[ii]]$lat <-
          model.grid$lat[clippings.subreg.list[[ii]]$offset.lon.lat[1]:
                         (clippings.subreg.list[[ii]]$offset.lon.lat[1] +
                          clippings.subreg.list[[ii]]$count.lon.lat[1] - 1),
                         clippings.subreg.list[[ii]]$offset.lon.lat[2]:
                         (clippings.subreg.list[[ii]]$offset.lon.lat[2] +
                          clippings.subreg.list[[ii]]$count.lon.lat[2] - 1)]
      }

      ## NAME SUBREGIONS
      ##  In case subregions already have a name, concatenate it with subregion
      ##  name, else just give single subregion name.
      has.names <- !is.null(names(clippings.subreg.list))
      if (has.names){
        names(clippings.subreg.list) <- paste(subregion.name,
                                              names(clippings.subreg.list),
                                              sep=".")
      } else {
        names(clippings.subreg.list) <- subregion.name
      }

    }

    output.list <- c(output.list, clippings.subreg.list)
  }

  return(output.list)
}



GetSubregionPolygons <- function(subregion, subregion.type, ...) {
  ## Get polygons (class "sp") out of a shapefile or rectangle.
  ##
  ## Args:
  ##   subregion: Single subregion from user.input
  ## subregion.type: Character indicating type of subregion
  ##                 currently:
  ##                  "rectangle" and "shapefile"
  ##                  (one could also calculate the convex hull of subregions as
  ##                   done previously...)
  ##
  ## Returns:
  ##   List of SpatialPolygons embedding subregions.
  ##
  ## History:
  ##   2010-10-25 | Original code (thm)
  ##   2010-11-16 | Change 'GetSubregionIndices' to 'GetSubregionPolygons' (thm)
  ##   2010-11-22 | Change shapefile and single.subregion user.input tags (thm)
  ##   2010-11-24 | Reading rectangles and shapefiles only. Parts of function
  ##                transferred to GetSubregionShapes (thm)
  ##   2012-08-29 | introducing: "category.only" to limit the evaluated
  ##                subregions to a given number (msu)

  ## If rectangle.
  ## Western, Eastern, Northern and Southern extensions
  if (subregion.type == "rectangle") {
    ## rectangle
   poly.sp <- Rectangle2Polygon(subregion)
  }
  ## If shapefile.
  if (subregion.type == "shapefile") {

    ## shapefile
    dirname <- subregion$dirname
    filename <- subregion$filename
    subregion.names <- subregion$subregionnames
    projection <- subregion$projection
    category <- subregion$category.variable
    category.label <- subregion$category.label
    category.only <- subregion$category.only
    poly.sp <- ShapeFile2Polygon(filename, dirname, proj = projection,
                                 category = category,
                                 category.label = category.label,
                                 category.only = category.only)

    ## if declared in user.input, specific subregionnames are given
    if (!is.null(subregion.names))
      names(poly.sp) <- subregion.names
  }

  return(poly.sp)
}



Rectangle2Polygon <- function(rec) {
  ## Makes a polygon (class "sp") out of western, eastern, northern and
  ## southern extansion rectangle. We need package "sp" for this function.
  ##
  ## Args:
  ##   rec: Western, Eastern, Northern and Southern extensions.
  ##        (vector of 'lon, lon, lat, lat')
  ##
  ## Returns:
  ##   List of single class SpatialPolygons rectangle.
  ##
  ## History:
  ##   2010-10-25 | Original code (thm)
  ##   2010-11-17 | Return value is list with single output (thm)

  ## generating polygon from input
  lon.poly <- c(rec[1], rec[1],rec[2],rec[2],rec[1])
  lat.poly <- c(rec[3], rec[4],rec[4],rec[3],rec[3])
  ## polygon sp class
  poly <- Polygons(list(Polygon(cbind(lon.poly, lat.poly))), "poly")
  poly.sp <- SpatialPolygons(list(poly))
  poly.sp <- list(poly.sp)

  return(poly.sp)
}



ShapeFile2Polygon <- function(file.shp, dir.shp, proj = NULL, category = NULL,
                              category.label = NULL, category.only = NULL) {
  ## Makes a polygon (class "sp") out of a shapefile.
  ## We need packages "rgdal" and "sp" for this function.
  ##
  ## Args:
  ##   file.shp: Filename of shapefile, without ending.
  ##   dir.shp: Directory of shapefile
  ##   category: category variable in shapefile. Polygons in one categoriy
  ##              will be sticked together.
  ##   category.label: Named vector for categories.
  ##
  ## Returns:
  ##   List of class SpatialPolygons polygons.
  ##
  ## History:
  ##   2010-10-25 | Original code (thm)
  ##   2010-11-16 | read multiple "Polygons" from single shapefile (thm)
  ##   2011-01-13 | case we have "SpatialLinesDataFrame" from readOGR as
  ##                shapefile...
  ##   2011-04-21 | set shapefile proj4string to WGS84 manualy if not set
  ##                (geh, thm)
  ##   2011-06-01 | optional category and category.label keyword (thm)
  ##   2011-06-29 | bugfix extracting multiple "Polygon" objects from within
  ##                "Polygons" (thm)
  ##   2012-08-29 | introducing: "category.only" to limit the evaluated
  ##                subregions to a given number (msu)

  ## generating sp polygon from input shapefile
  cat("    reading in shapefile...\n")
  shp.sp <- rgdal::readOGR(dir.shp, file.shp, verbose = FALSE)

  ## error if proj4string not set
  if (is.null(proj)) {
    if (is.na(proj4string(shp.sp)))
      stop("NO PROJECTION SPECIFIED IN SHAPEFILE, PLEASE CONTACT THE SHAPEFILE PROVIDER OR SET A PROJECTION IN user.input LIKE projection = \"+proj=longlat +ellps=WGS84\"")
  } else {
    proj4string(shp.sp) <- CRS(proj)
  }
  ## transform to WGS84 lat and lon
  llCRS <- CRS("+init=epsg:4326")
  shp.sp <- spTransform(shp.sp, llCRS)

### readOGR can return a SpatialPolygonsDataFrame object or a
### SpatialLinesDataFrame object...
### extract "Polygons" slot from "SpatialPolygons" object

  ## TODO, if categories given by user.input.... else read list element by
  ## list element as earlier (NA, so wie jetzt aber categorien mit rownames
  ## ersetzen)
    if (!is.null(category)) {
      categ <- shp.sp@data[[category]]
    } else {
      categ <- as.numeric(row.names(shp.sp@data))
    }
  categ.levels <- unique(sort(categ))

  #### msu HACK ####
  ## init output list
  if ( !is.null(category.only) ) {
    polygons <- vector("list", length(category.only))
  } else {
    polygons <- vector("list", length(categ.levels))
  }
  #### /msu HACK ####

  ## case of SpatialPolygonsDataFrame
  if (class(shp.sp) == "SpatialLinesDataFrame") {

    if (length(shp.sp@lines) != length(categ))
      stop("NUMBER OF CATEGORIES WITHIN SHAPEFILE DOESNT MATCH NUMBER OF ENTRIES")
    ## init subreg names to be defined later
    subregion.names <- NULL

    ## convert "SpatialLinesDataFrame" to list of "polygons"
    lines <- sapply(slot(shp.sp, "lines"), function(x) x)
    ## concatenate polygons within the same catogory to a 'Polygons' object
    for (ii in seq(along = categ.levels)) {
      categ.value <- categ.levels[ii]
      which.polys.in.categ <- which(categ %in% categ.value)
      ## set name for catogory
      if (!is.null(category.label) > 0)
        ID.name <- category.label[categ.value == names(category.label)]
      else if (is.null(category.label) > 0) ## && is.factor(categ.levels) (old?)
        ID.name <-  categ.value
      else
        ID.name = ii - 1

      ## define helperfunctions to extract lon/lat coordinates from a "Lines"
      ## object and to create a list of Polygon objects from
      ## matrix of lon/lat coordinates
      Lines2coordmatrix <- function(x){
        ## Gets list of lon-lat matrices out of a Lines object (containng
        ## (a list of) Line objects).
        ## x: Lines object (package sp)
        lines.list <- slot(x, "Lines")
        y <- lapply(lines.list, function(z) slot(z, "coords"))
        return(y)
      }
      coordmatrix2Polygon <- function(x) {
        ## generates list of Polygon objects out of list of coordinate matrices
        ## x: list of lon lat matrices
        y <- sapply(x, Polygon)
        return(y)
      }

      ## converts list of 'Lines' sp-objects to a list of 'Polygon' objects
      polys.in.categ <- sapply(which.polys.in.categ, function(x)
                              coordmatrix2Polygon(Lines2coordmatrix(lines[[x]]))
                               )
      ## converts this list of 'Polygon' objects to a single 'Polygons' object
      polygons[[ii]] <- Polygons(polys.in.categ, ID = ID.name)
      ## Names of current subregion (string or just ID, depends on categories)
      subregion.names <- c(subregion.names, as.character(ID.name))
    }
  }
  ## case of SpatialPolygonsDataFrame
  if (class(shp.sp) == "SpatialPolygonsDataFrame") {

    ## init subreg names to be defined later
    subregion.names <- NULL

    #### msu HACK ####
    if ( !is.null(category.only) ) {
      categ.levels <- category.only
    }
    #### /msu HACK ####

    for (ii in seq(along = categ.levels)) {
      id <- ii - 1
      categ.value <- categ.levels[ii]
      which.polys.in.categ <- which(categ %in% categ.value)
      ## set name for catogory
      if (!is.null(category.label) > 0)
        ID.name <- category.label[categ.value == names(category.label)]
      else if (is.null(category.label) > 0) ## && is.factor(categ.levels) (old?)
        ID.name <- categ.value
      else
        ID.name = id

      ## extracts polygons belonging to the same class to a list
      polys.in.categ <- unlist(lapply(which.polys.in.categ, function(x)
                                      slot(slot(shp.sp, "polygons")[[x]],
                                           "Polygons")
                                      )
                             )

      ## converts this list of 'Polygon' objects to a single 'Polygons' object
      polygons[[ii]] <- Polygons(polys.in.categ, ID = id)
      ## Names of current subregion (string or just ID, depends on categories)
      subregion.names <- c(subregion.names, as.character(ID.name))
    }
  }
  ##  case of SpatialPointsDataFrame (didnt occur until now)
  if (class(shp.sp) == "SpatialPointsDataFrame") {
    stop("SpatialPointsDataFrame OBJECT NOT YET IMPLEMENTED IN ShapeFile2Polygon")
  }

  ## get the "ID" of each "Polygons"
  ids <- sapply(polygons, function(x) slot(x, "ID"))
  ## generate list of "SpatialPolygons" for each "Polygons"
  sp.list <- lapply(polygons, function(x) SpatialPolygons(list(x),
                                                          proj4string=llCRS))
  names(sp.list) <- subregion.names
  cat("   ", length(ids), "polygons created by shapefile\n")

  return(sp.list)
}



SubregionFile2ModelGrid <- function(subregion, target.grid) {
  ## Interpolates subregionsfile to climate model grid via nearest neighbor.
  ## First, the subregion will get an "envelope" with NAs, so nearest neighbor
  ## gets the right region (see BuildSubregionEnvelope).
  ## Needs package "class" to run nearest neighbor.
  ##
  ## Args:
  ##   subregion: subregion tag from user.input containing subregionsfile info.
  ##   target.grid: Full model grid file name (with directory).
  ##
  ## Returns:
  ##   List of mask-matrices on model grid with NAs indicating region lying
  ##   outside the specified subregion.
  ##
  ## History:
  ##   2010-11-24 | Original code (thm)

  ## extract information from subregion-list (function argument)
  subreg.grid.file <- paste(subregion$grid.dir, subregion$grid.file, sep="/")
  subreg.file <- paste(subregion$subreg.dir, subregion$subreg.file, sep="/")
  mask.name <- subregion$mask.name
  mask.values <- subregion$mask.value

  ## set default mask variable name in nc file
  if(is.null(mask.name))
    mask.name <- "mask"

  ## read in lat and lon from subregion
  cat("  ")
  subreg.grid <- ReadGridFile(subreg.grid.file)
  subreg.lon <- subreg.grid$lon
  subreg.lat <- subreg.grid$lat

  ## get lat and lon from target grid (climate model)
  target.xy <-  matrix(c(c(target.grid$lon), c(target.grid$lat)), ncol=2)
  target.dims <- target.grid$dims

  ## initialize list to be returned
  model.subregions.list <- vector("list", length(mask.values))

  ## read subregion after subregion
  for (mask.count in seq(along = mask.values)) {
    ## get subregion mask value
    single.mask.variable <- mask.values[mask.count]
    ## read in subregion mask indices
    subreg <- ReadGridFile(subreg.file, mask.name)
    subreg <- subreg$mask
    is.na(subreg) <- which(subreg != single.mask.variable)
    ## build envelope so that subregion is surrounded by NAs, else nearest
    ## neighbor calculates the wrong subregions on model grid
    ## eg see:
    ##   stop()
    ##   image(subreg)
    ##   subreg.fac <- c(subreg)
    ##   subreg.fac[is.na(subreg.fac)] <- "outside"
    ##   xy <-  matrix(c(c(subreg.lon), c(subreg.lat)), ncol=2)
    ##   result <- knn1(xy, target.xy, subreg.fac)
    ##   is.na(result) <- as.numeric(which(result == "outside"))
    ##   image(matrix(as.numeric(result), ncol=target.dims[["lat"]]))
    grid.env <- BuildSubregionEnvelope(subreg = subreg,
                                       lon = subreg.lon, lat = subreg.lat)
    subreg.env <- grid.env$subreg
    lon.env <- grid.env$lon
    lat.env <- grid.env$lat
    xy.env <-  matrix(c(c(lon.env), c(lat.env)), ncol=2)
    subreg.env.fac <- c(subreg.env)
    ## NAs have to be flaged as "outside",as nearest neighbor does not allow NAs
    subreg.env.fac[is.na(subreg.env.fac)] <- "outside"

    ## nearest neighbor with enveloped grid and subregion
    ## now the subregion data will be on the specific model grid
    cat("    nearest neighbor to model grid for subregion mask",
        single.mask.variable, "\n")
    result.env <- class::knn1(xy.env, target.xy, subreg.env.fac)
    ## get back "outside" to NAs and reform results to matrix
    is.na(result.env) <- which(result.env == "outside")
    result.env <- matrix(as.numeric(result.env), ncol = target.dims[["lat"]])

    ## put subregion on model grid to list to be returned
    model.subregions.list[[mask.count]] <- result.env
   }

### subregions namings

  ## get indicator variables
  has.more.than.one.subregion <- length(mask.values) > 1
  has.subregionnames <- (!is.null(subregion$subregionnames))

  ## in case we have more than one subregion, we have to name them in
  ## order to know which subregion is which (take the mask.values as names)
  if (has.more.than.one.subregion)
    names(model.subregions.list) <- mask.values
  ## in case we have subregionnames in user.input file, name after them
  if (has.subregionnames)
    names(model.subregions.list) <- subregion$subregionnames

  return(model.subregions.list)
}



GetIndicesFromPolygon <- function(poly.sp, coordinates, wrap.to = "180",
                                  area.fraction = FALSE) {
  ## Get offset and count for reading in rectangular array from NetCDF
  ## file and index mask for clipping polygon 'poly.sp' from this rectangle.
  ##
  ## Args:
  ##   poly.sp: SpatialPolygon class object ("sp" package)
  ##   coordinates: List of lat lon matrices obtained from "ReadGridFile"
  ##   wrap.to: workaround for -180 to 180 lon-problem
  ##   area.fraction: boolean, if true masking of the cells is done by the
  ##      area fraction of cell and overlaying polygon
  ##
  ## Returns:
  ##   List containing 'offset', 'count' for subregion, 'mask' for exact
  ##   clipping of polygon.
  ##
  ## History:
  ##   2010-10-25 | Original code (thm)
  ##   2010-11-17 | Return value is list with single output (thm)
  ##   2011-02-15 | wrap to 180 and to 360 added (thm)
  ##   2011-02-16 | wrap to 180 and to 360 outsourced as own function (thm)
  ##   2011-07-08 | masking cells with area fraction of overlaying polygon
  ##                (mir & thm)
  ##   2014-10-10 | changed deprecated function "overlay" to "over" (thm)

  nc.dims <- coordinates$dims
  ## generate spatial object model grid
  lon <- c(coordinates$lon)
  lat <- c(coordinates$lat)
  lonlat.proj <- "+proj=longlat +ellps=WGS84"
  whole.region.sp <- SpatialPoints(cbind(lon, lat),
                                   proj4string = CRS(lonlat.proj))

  proj4string(poly.sp) <- CRS(lonlat.proj)
  in.poly <- sp::over(whole.region.sp, poly.sp)

  ## when all NA
  if (all(is.na(in.poly)))
    stop("NO GRIDFILE COORDINATES LIE WITHIN SPECIFIED POLYGON.")

  ## in.poly must be a vector. if we use SpatialPolygonsDataFrame objects,
  ## 'overlay' generates a data.frame instead of a simple vector
  if (is.data.frame(in.poly))
    in.poly <- in.poly$ID

  ## get to matrix form again
  in.poly <- matrix(in.poly, ncol = nc.dims[["lat"]],
                    nrow = nc.dims[["lon"]])

  ## Retrieve spatial count and offset
  cat("  ")
  ret <- GetOffsetCountAndMask(in.poly, nc.dims)

  return(ret)
}



GetOffsetCountAndMask <- function(mask.mat, dims) {
  ## Get offset and count for reading in rectangular array from NetCDF
  ## file and index mask for clipping subregion from this rectangule.
  ##
  ## Args:
  ##   mask.mat: Matrix with NAs indicating region lying outside subregion.
  ##   dims: Original dim size of netcdf file.
  ##         In fact its not really necessary, because dims == dim(mask.mat)?
  ##
  ## Returns:
  ##   List containing 'offset', 'count' for subregion, 'mask' for exact
  ##   clipping of subregion.
  ##
  ## History:
  ##   2010-10-25 | Original code (thm)

  cat("  ")

  ## get the indices of subregion
  mask.inside.index  <- which(!is.na(mask.mat))

  ## look for max and min index of lon and lat values in subregion.
  ## now not array anymore, but coerced by 'which' to vector in R,
  ## therefore one has to calculate the 2 dimensional indices out of
  ## 1-dimensional vector indices.
  ## The lat value (y-axis) index is obtained by dividing by length of lon.
  ## The lon value (x-axis) by ounting the "rest" of this former division
  ## with a modulo operation.
  latitude.index.in.subregion <- ceiling(mask.inside.index / dims["lon"])
  longitude.index.in.subregion <- mask.inside.index %% dims["lon"]

  ## Longitude: case modulo = 0 means LAST row, i.e. add dims["lon"]
  ## (in general we have dims["lon"] rows and dims["lat"] columns)
  ## To get column index, we divide by dims["lon"].
  longitude.index.in.subregion[longitude.index.in.subregion == 0] <- dims["lon"]

  ## get the min-max index boundaries
  min.longitude.index.in.subregion <- min(longitude.index.in.subregion)
  max.longitude.index.in.subregion <- max(longitude.index.in.subregion)
  min.latitude.index.in.subregion <- min(latitude.index.in.subregion)
  max.latitude.index.in.subregion <- max(latitude.index.in.subregion)

  ## calculate offset and count
  offset <- c(min.longitude.index.in.subregion, min.latitude.index.in.subregion)
  count <- c(max.longitude.index.in.subregion, max.latitude.index.in.subregion) - offset + 1

  ## get mask array for setting to NA
  mask.subregion.mat <- mask.mat[offset[1]:(offset[1] - 1 + count[1]),
                                 offset[2]:(offset[2] - 1 + count[2])]

  ## indizes for pixels lying outside the region of interest
  mask.outside.index <- which(is.na(mask.subregion.mat))

### return list with offset, count and index mask for values outside subregion
### within offset-count region
  ret <- list(offset.lon.lat = offset,
              count.lon.lat = count,
              mask.indices.outside.subregion = mask.outside.index)

  return(ret)
}



BuildSubregionEnvelope <- function(subreg, lon, lat) {
  ## Creates a NA-column or NA-row around subregion matrix if at least one value
  ## on the border has not an NA value. This is done as a step before using
  ## nearest neighbor, to obtain a matrix with first/last row/column filled with
  ## NAs to close the subregion (else nearest neighbor gets wrong values
  ## - see example in SubregionFile2ModelGrid code).
  ## lon and lat matrices are being extended by corresponding rows/columns,
  ## by taking the diffrence between the neighboring columns or rows (that is a
  ## bit illegal for irregular grids).
  ##
  ## Args:
  ##   subreg: Matrix containing NAs indicating values outside the subregion.
  ##   lon: Matrix containing subregion longitudes.
  ##   lat: Matrix containing subregion latitudes.
  ##
  ## Returns:
  ##   List with extended subreg, lon and lat matrices.
  ##
  ## History:
  ##   2010-11-24 | Original code (thm)
  ##
  ## TODO:
  ##   Check if regular grid. if not, give a message: "envelope only
  ##   approxmation"

  ## error if dimensions do not match
  if (!all(dim(subreg) == c(dim(lat), dim(lon))))
    stop("SUBREGION DIMENSION DON'T MATCH GRID DIMENSIONS.")
  ## all input arguments must be a matrix
  if (!is.matrix(subreg))
    stop("ARGUMENTS MUST BE MATRICES.")

  ## get last matrix row/column index
  last.row <- dim(subreg)[1]
  last.col <- dim(subreg)[2]

  ## declare envelope subregions lon and lat
  subreg.fake.mat <- subreg
  lon.fake.mat <- lon
  lat.fake.mat <- lat

  ## northern borderline (i.e. last column in R)
  ## if there is at least one element != NA then make envelope
  if (any(!is.na(subreg.fake.mat[, last.col]))) {
    ## add 0 to subregionboarder
    subreg.fake.mat <- cbind(subreg.fake.mat, NA)
    ## calculate latitude and longitude extension by taking the difference
    fake.north.lat <- 2*lat.fake.mat[, last.col] - lat.fake.mat[, last.col - 1]
    fake.north.lon <- 2*lon.fake.mat[, last.col] - lon.fake.mat[, last.col - 1]
    ## add fake northern lat to lat matrix
    lat.fake.mat <-  cbind(lat.fake.mat, fake.north.lat)
    ## add fake northern lon to lon matrix
    lon.fake.mat <-  cbind(lon.fake.mat, fake.north.lon)
  }
  ## southern borderline (i.e. first column in R)
  ## if there is at least one element != NA then make envelope
  if (any(!is.na(subreg.fake.mat[, 1]))) {
    subreg.fake.mat <- cbind(NA, subreg.fake.mat)
    fake.south.lat <- 2*lat.fake.mat[, 1] - lat.fake.mat[, 2]
    fake.south.lon <- 2*lon.fake.mat[, 1] - lon.fake.mat[, 2]
    ## add fake southern lat to lat matrix
    lat.fake.mat <-  cbind(lat.fake.mat, fake.south.lat)
    ## add fake southern lon to lon matrix
    lon.fake.mat <-  cbind(lon.fake.mat, fake.south.lon)
  }
  ## eastern boarderline (i.e. last row in R)
  ## if there is at least one element != NA then make envelope
  if (any(!is.na(subreg.fake.mat[last.row, ]))){
    subreg.fake.mat <- rbind(subreg.fake.mat, NA)
    fake.east.lat <- 2*lat.fake.mat[last.row, ] - lat.fake.mat[last.row - 1, ]
    fake.east.lon <- 2*lon.fake.mat[last.row, ] - lon.fake.mat[last.row - 1, ]
    ## add fake eastern lat to lat matrix
    lat.fake.mat <-  rbind(lat.fake.mat, fake.east.lat)
    ## add fake eastern lon to lon matrix
    lon.fake.mat <-  rbind(lon.fake.mat, fake.east.lon)
  }
  ## western boarderline (i.e. first row in R)
  ## if there is at least one element != NA then make envelope
  if (any(!is.na(subreg.fake.mat[1, ]))) {
    subreg.fake.mat <- rbind(NA, subreg.fake.mat)
    fake.west.lat <- 2*lat.fake.mat[1, ] - lat.fake.mat[2, ]
    fake.west.lon <- 2*lon.fake.mat[1, ] - lon.fake.mat[2, ]
    ## add fake western lat to lat matrix
    lat.fake.mat <-  rbind(fake.west.lat, lat.fake.mat)
    ## add fake western lon to lon matrix
    lon.fake.mat <-  rbind(fake.west.lon, lon.fake.mat)
  }
  ## erase dimnames caused by rbind (plastical surgery)
  dimnames(lon.fake.mat) <- NULL
  dimnames(lat.fake.mat) <- NULL
  ## define return value
  subreg.lon.lat <- list(subreg = subreg.fake.mat,
                         lon = lon.fake.mat,
                         lat = lat.fake.mat)
  return(subreg.lon.lat)
}



BoundingPoly <- function(my.grid, i, j) {
  ## get the cell bounds to index i/j of the given grid
  ## as spatial polygon
  ## if the cell lies inside the grid, the bounding box
  ## can easily be calculated using the neighbor cells
  ## cell at the grid border are computed by translating
  ## the bounding box of the nearest in-grid neighbor cell
  ## to the grid boundary
  ##
  ## Args:
  ##   my.grid: dataframe containing lon/lat matrix
  ##   i: integer row index
  ##   j: integer column index
  ##
  ## Returns:
  ##   SpatialPolygon containing cell boundary
  ##
  ## History:
  ##   2011-07-08 | Original code (mir & thm)
  ##   2011-10-17 | introduced correct special cases
  ##              | for corners and boundaries (geh)


  ## get grid dimensions
  dims <- dim(my.grid$lon)

  ## error if i or j are outside grid
  if( (i > dims[1]) || (j > dims[2]) ) {
    cat("ERROR in BoundingPoly: GRID INDICES ARE OUTSIDE THE GRID \n")
  }

  ## inside the grid -- calc the corners using the eight neighbours
  ii <- max(2, i)
  ii <- min(dims[1]-1, ii)
  jj <- max(2, j)
  jj <- min(dims[2]-1, jj)

  sw.lon <- my.grid$lon[ii-1, jj-1]
  sw.lat <- my.grid$lat[ii-1, jj-1]
  s.lon <- my.grid$lon[ii, jj-1]
  s.lat <- my.grid$lat[ii, jj-1]
  se.lon <- my.grid$lon[ii+1, jj-1]
  se.lat <- my.grid$lat[ii+1, jj-1]

  w.lon <- my.grid$lon[ii-1, jj]
  w.lat <- my.grid$lat[ii-1, jj]
  c.lon <- my.grid$lon[ii, jj]
  c.lat <- my.grid$lat[ii, jj]
  e.lon <- my.grid$lon[ii+1, jj]
  e.lat <- my.grid$lat[ii+1, jj]

  nw.lon <- my.grid$lon[ii-1, jj+1]
  nw.lat <- my.grid$lat[ii-1, jj+1]
  n.lon <- my.grid$lon[ii, jj+1]
  n.lat <- my.grid$lat[ii, jj+1]
  ne.lon <- my.grid$lon[ii+1, jj+1]
  ne.lat <- my.grid$lat[ii+1, jj+1]

  ## special case for south-eastern corner
  if ( (i == 1 && j == 1) ) {
    trans.lon <- my.grid$lon[i, j] - c.lon
    trans.lat <- my.grid$lat[i, j] - c.lat
    n.lon <- w.lon
    n.lat <- w.lat
    ne.lon <- c.lon
    ne.lat <- c.lat
    e.lon <- s.lon
    e.lat <- s.lat
    c.lon <- my.grid$lon[i, j]
    c.lat <- my.grid$lat[i, j]

    sw.lon <- sw.lon+trans.lon
    sw.lat <- sw.lat+trans.lat
    s.lon <- s.lon+trans.lon
    s.lat <- s.lat+trans.lat
    se.lon <- se.lon+trans.lon
    se.lat <- se.lat+trans.lat

    w.lon <- w.lon+trans.lon
    w.lat <- w.lat+trans.lat

    nw.lon <- nw.lon+trans.lon
    nw.lat <- nw.lat+trans.lat
  }

  ## special case for north-eastern corner
  if ( (i == 1 && j == dims[2]) ) {
    trans.lon <- my.grid$lon[i, j] - c.lon
    trans.lat <- my.grid$lat[i, j] - c.lat
    s.lon <- w.lon
    s.lat <- w.lat
    se.lon <- c.lon
    se.lat <- c.lat
    e.lon <- n.lon
    e.lat <- n.lat
    c.lon <- my.grid$lon[i, j]
    c.lat <- my.grid$lat[i, j]

    sw.lon <- sw.lon+trans.lon
    sw.lat <- sw.lat+trans.lat

    w.lon <- w.lon+trans.lon
    w.lat <- w.lat+trans.lat

    nw.lon <- nw.lon+trans.lon
    nw.lat <- nw.lat+trans.lat
    n.lon <- n.lon+trans.lon
    n.lat <- n.lat+trans.lat
    ne.lon <- ne.lon+trans.lon
    ne.lat <- ne.lat+trans.lat
  }

  ## special case for north-western corner
  if ( (i == dims[1] && j == dims[2]) ) {
    trans.lon <- my.grid$lon[i, j] - c.lon
    trans.lat <- my.grid$lat[i, j] - c.lat
    s.lon <- e.lon
    s.lat <- e.lat
    sw.lon <- c.lon
    sw.lat <- c.lat
    w.lon <- n.lon
    w.lat <- n.lat
    c.lon <- my.grid$lon[i, j]
    c.lat <- my.grid$lat[i, j]

    se.lon <- se.lon+trans.lon
    se.lat <- se.lat+trans.lat

    e.lon <- e.lon+trans.lon
    e.lat <- e.lat+trans.lat

    nw.lon <- nw.lon+trans.lon
    nw.lat <- nw.lat+trans.lat
    n.lon <- n.lon+trans.lon
    n.lat <- n.lat+trans.lat
    ne.lon <- ne.lon+trans.lon
    ne.lat <- ne.lat+trans.lat
  }

  ## special case for south-western corner
  if ( (i == dims[1] && j == 1) ) {
    trans.lon <- my.grid$lon[i, j] - c.lon
    trans.lat <- my.grid$lat[i, j] - c.lat
    n.lon <- e.lon
    n.lat <- e.lat
    nw.lon <- c.lon
    nw.lat <- c.lat
    w.lon <- s.lon
    w.lat <- s.lat
    c.lon <- my.grid$lon[i, j]
    c.lat <- my.grid$lat[i, j]

    sw.lon <- sw.lon+trans.lon
    sw.lat <- sw.lat+trans.lat
    s.lon <- s.lon+trans.lon
    s.lat <- s.lat+trans.lat
    se.lon <- se.lon+trans.lon
    se.lat <- se.lat+trans.lat

    e.lon <- e.lon+trans.lon
    e.lat <- e.lat+trans.lat

    ne.lon <- ne.lon+trans.lon
    ne.lat <- ne.lat+trans.lat
  }

  ## special case for western boundary
  if ( (i == 1) && (j >= 2) && (j <= dims[2]-1) ) {
    sw.lon <- my.grid$lon[i, j-1] +
      my.grid$lon[i, j-1] - my.grid$lon[i+1, j-1]
    sw.lat <- my.grid$lat[i, j-1] +
       my.grid$lat[i, j-1] - my.grid$lat[i+1, j-1]
    s.lon <- my.grid$lon[i, j-1]
    s.lat <- my.grid$lat[i, j-1]
    se.lon <- my.grid$lon[i+1, j-1]
    se.lat <- my.grid$lat[i+1, j-1]

    w.lon <- my.grid$lon[i, j] +
      my.grid$lon[i, j] - my.grid$lon[i+1, j]
    w.lat <- my.grid$lat[i, j] +
      my.grid$lat[i, j] - my.grid$lat[i+1, j]
    c.lon <- my.grid$lon[i, j]
    c.lat <- my.grid$lat[i, j]
    e.lon <- my.grid$lon[i+1, j]
    e.lat <- my.grid$lat[i+1, j]

    nw.lon <- my.grid$lon[i, j+1] +
      my.grid$lon[i, j+1] - my.grid$lon[i+1, j+1]
    nw.lat <-  my.grid$lat[i, j+1] +
      my.grid$lat[i, j+1] - my.grid$lat[i+1, j+1]
    n.lon <- my.grid$lon[i, j+1]
    n.lat <- my.grid$lat[i, j+1]
    ne.lon <- my.grid$lon[i+1, j+1]
    ne.lat <- my.grid$lat[i+1, j+1]
  }

  ## special case for northern boundary
  if ( (j == dims[2]) && (i >= 2) && (i <= dims[1]-1) ) {
    sw.lon <- my.grid$lon[i-1, j-1]
    sw.lat <- my.grid$lat[i-1, j-1]
    s.lon <- my.grid$lon[i, j-1]
    s.lat <- my.grid$lat[i, j-1]
    se.lon <- my.grid$lon[i+1, j-1]
    se.lat <- my.grid$lat[i+1, j-1]

    w.lon <- my.grid$lon[i-1, j]
    w.lat <- my.grid$lat[i-1, j]
    c.lon <- my.grid$lon[i, j]
    c.lat <- my.grid$lat[i, j]
    e.lon <- my.grid$lon[i+1, j]
    e.lat <- my.grid$lat[i+1, j]

    nw.lon <- my.grid$lon[i-1, j] +
      my.grid$lon[i-1, j] - my.grid$lon[i-1, j-1]
    nw.lat <-  my.grid$lat[i-1, j] +
      my.grid$lat[i-1, j] - my.grid$lat[i-1, j-1]
    n.lon <- my.grid$lon[i, j] +
      my.grid$lon[i, j] - my.grid$lon[i, j-1]
    n.lat <- my.grid$lat[i, j] +
      my.grid$lat[i, j] - my.grid$lat[i, j-1]
    ne.lon <- my.grid$lon[i+1, j] +
      my.grid$lon[i+1, j] - my.grid$lon[i+1, j-1]
    ne.lat <- my.grid$lat[i+1, j] +
      my.grid$lat[i+1, j] - my.grid$lat[i+1, j-1]
  }

  ## special case for eastern boundary
  if ( (i == dims[1]) && (j >= 2) && (j <= dims[2]-1) ) {
    sw.lon <- my.grid$lon[i-1, j-1]
    sw.lat <- my.grid$lat[i-1, j-1]
    s.lon <- my.grid$lon[i, j-1]
    s.lat <- my.grid$lat[i, j-1]
    se.lon <- my.grid$lon[i, j-1] +
      my.grid$lon[i, j-1] - my.grid$lon[i-1, j-1]
    se.lat <- my.grid$lat[i, j-1] +
       my.grid$lat[i, j-1] - my.grid$lat[i-1, j-1]

    w.lon <- my.grid$lon[i-1, j]
    w.lat <- my.grid$lat[i-1, j]
    c.lon <- my.grid$lon[i, j]
    c.lat <- my.grid$lat[i, j]
    e.lon <- my.grid$lon[i, j] +
      my.grid$lon[i, j] - my.grid$lon[i-1, j]
    e.lat <- my.grid$lat[i, j] +
      my.grid$lat[i, j] - my.grid$lat[i-1, j]

    nw.lon <- my.grid$lon[i-1, j+1]
    nw.lat <-  my.grid$lat[i-1, j+1]
    n.lon <- my.grid$lon[i, j+1]
    n.lat <- my.grid$lat[i, j+1]
    ne.lon <- my.grid$lon[i, j+1] +
      my.grid$lon[i, j+1] - my.grid$lon[i-1, j+1]
    ne.lat <- my.grid$lat[i, j+1] +
      my.grid$lat[i, j+1] - my.grid$lat[i-1, j+1]
  }

  ## special case for southern boundary
  if ( (j == 1) && (i >= 2) && (i <= dims[1]-1) ) {
    sw.lon <- my.grid$lon[i-1, j] +
      my.grid$lon[i-1, j] -  my.grid$lon[i-1, j+1]
    sw.lat <- my.grid$lat[i-1, j] +
      my.grid$lat[i-1, j] -  my.grid$lat[i-1, j+1]
    s.lon <- my.grid$lon[i, j] +
      my.grid$lon[i, j] -  my.grid$lon[i, j+1]
    s.lat <- my.grid$lat[i, j] +
      my.grid$lat[i, j] -  my.grid$lat[i, j+1]
    se.lon <- my.grid$lon[i+1, j] +
      my.grid$lon[i+1, j] -  my.grid$lon[i+1, j+1]
    se.lat <- my.grid$lat[i+1, j] +
      my.grid$lat[i+1, j] -  my.grid$lat[i+1, j+1]

    w.lon <- my.grid$lon[i-1, j]
    w.lat <- my.grid$lat[i-1, j]
    c.lon <- my.grid$lon[i, j]
    c.lat <- my.grid$lat[i, j]
    e.lon <- my.grid$lon[i+1, j]
    e.lat <- my.grid$lat[i+1, j]

    nw.lon <- my.grid$lon[i-1, j+1]
    nw.lat <- my.grid$lat[i-1, j+1]
    n.lon <- my.grid$lon[i, j+1]
    n.lat <-  my.grid$lat[i, j+1]
    ne.lon <-  my.grid$lon[i+1, j+1]
    ne.lat <- my.grid$lat[i+1, j+1]
  }

  ## compute the bounding box by caluclating the midpoint of
  ## the corresponding neighbor-centroids
  sw.point <- Polygon(cbind(c(sw.lon, s.lon, c.lon, w.lon, sw.lon),
                            c(sw.lat, s.lat, c.lat, w.lat, sw.lat)))@labpt
  se.point <- Polygon(cbind(c(se.lon, e.lon, c.lon, s.lon, se.lon),
                            c(se.lat, e.lat, c.lat, s.lat, se.lat)))@labpt
  ne.point <- Polygon(cbind(c(ne.lon, n.lon, c.lon, e.lon, ne.lon),
                            c(ne.lat, n.lat, c.lat, e.lat, ne.lat)))@labpt
  nw.point <- Polygon(cbind(c(nw.lon, w.lon, c.lon, n.lon, nw.lon),
                            c(nw.lat, w.lat, c.lat, n.lat, nw.lat)))@labpt

  ## build the spatial polygon with uniq id (running index of the cell in c-notation)
  ret <-
    SpatialPolygons(list(Polygons(list(Polygon(cbind(c(sw.point[1], se.point[1],
                                                       ne.point[1],
                                                       nw.point[1],
                                                       sw.point[1]),
                                                     c(sw.point[2], se.point[2],
                                                       ne.point[2],
                                                       nw.point[2],
                                                       sw.point[2])))),
                                       (jj-1)*dims[1]+ii-1)))

  return (ret)
}



GetWeightingMaskFromPolygon <- function(poly.sp, grid.list) {
  ## overlay cells of grid with the shape polygon and
  ## calculate the area fraction of this overlay for each cell
  ##
  ## Args:
  ##   poly.sp: SpatialPolygons representing one subregion
  ##   grid.list: list holding grid cell-centers
  ##
  ## Returns:
  ##   List containing 'offset', 'count' for subregion, 'mask' for exact
  ##   clipping of polygon, weight (area fraction) of cells.
  ##
  ## History:
  ##   2011-07-08 | Original code (mir & thm)
  ##   2012-08-29 | simply take nearest neighbour gridpoint if everything
  ##                else fails (msu)
  ##   2014-10-10 | changed deprecated function "overlay" to "over" (thm)
  ##
  cat ("      calculating areal fractions\n")

  ## sort grid according to monotonously ascending longitudes
  sort.indices.list <- sort(grid.list$lon[ ,1], index.return=TRUE)
  sort.grid.list <- grid.list
  sort.grid.list$lon <- grid.list$lon[sort.indices.list$ix, ]
  sort.grid.list$lat <- grid.list$lat[sort.indices.list$ix, ]
  dims <- dim(sort.grid.list$lon)

  ## initialize mask
  mask = matrix(0, nrow=dims[1], ncol=dims[2])

  ## just consider cells overlaying the bbox of shape + 1 row/col on each side
  rec <- c(t(bbox(poly.sp)))
  poly.bbox <- Rectangle2Polygon(rec)[[1]]
  whole.region.sp <- SpatialPoints(cbind(c(sort.grid.list$lon), c(sort.grid.list$lat)))

  ## clipping subregions
  in.poly <- matrix(sp::over(whole.region.sp, poly.bbox,),
                    ncol=dims[2], nrow=dims[1])
  
  if ( any(is.finite(in.poly)) ) {

    ## add one row / column to the rectangle of the grid that's covered by
    ## the bounding box of the considered shape
    ind.list <- which(is.finite(in.poly), arr.ind=T)
    mini <- max(min(ind.list[,1])-1, 1)
    maxi <- min(max(ind.list[,1])+1, dims[1])
    minj <- max(min(ind.list[,2])-1, 1)
    maxj <- min(max(ind.list[,2])+1, dims[2])

    ## double for loop over suspicious cells
    for ( ii in mini:maxi ) {
      for ( jj in minj:maxj ) {
        my.cell <- BoundingPoly(sort.grid.list, ii, jj)
        my.cell.area <- gArea(my.cell)
        ## for loop over polygons in shape (workaround to the rgeos package
        ## deficits (intersection can only be done between two single polygons)
        for ( p in 1:length(poly.sp@polygons[[1]]@Polygons) ) {
          ## factor, if polygon is a hole-polygon --> area fraction must be
          ## subtracted
          hole.factor <- poly.sp@polygons[[1]]@Polygons[[p]]@hole * -2 + 1
          my.small.poly <-
            SpatialPolygons(list(Polygons(list((poly.sp)@polygons[[1]]@Polygons[[p]]),ID=p)))
          if ( my.small.poly@polygons[[1]]@Polygons[[1]]@area > 10^-6 ) {
            intsec <- gIntersection(my.cell, my.small.poly)
            if ( is.null(intsec) ) area.fraction <- 0.
            else {
              area.fraction <- gArea(intsec)/my.cell.area
              mask[ii, jj] <- mask[ii, jj] + area.fraction * hole.factor
            }
          }
        }
      }
    }
  } else {
    ## no overlapping grid box has been found, hurray!
    ## so, search for the nearest grid point and take its value.
    cat ("      ATTENTION: No overlapping grid box with subregion has been found!\n",
         "      Need to take nearest grid point\n")

    euclidian.distance <- matrix(spDistsN1(whole.region.sp,
                                           poly.bbox@polygons[[1]]@labpt,
                                           longlat = FALSE),
                                 nrow = dims[1], ncol = dims[2])

    min.index <- which(euclidian.distance == min(euclidian.distance),
                       arr.ind = TRUE)
    mask[min.index[1],min.index[2]] <- 1
  }

  ## resort mask array according to original grid
  mask <- mask[match(grid.list$lon[ ,1], sort.grid.list$lon[ ,1]), ]

  mask[which(mask > 1)] <- 1
  mask[which(mask < 0)] <- 0
  mask[which(mask == 0)] <- NA
  names(dims) <- c('lon', 'lat')
  ret <- GetOffsetCountAndMask(mask, dims)
  offset <- ret$offset.lon.lat
  count <- ret$count.lon.lat
  ret$weight <- mask[offset[1]:(offset[1]+count[1]-1),
                     offset[2]:(offset[2]+count[2]-1)]

  return (ret)
}
