# Note that reg.finalizer does not finalize objects
# at the end of an R session. This could be a problem.
.setCollectorFun <- function(object, fun) {

  if (is.null(fun)) fun <- function(obj) obj
  reg.finalizer(object, fun, onexit=TRUE)

}

assertClass <- function(object, class) {
  
  if (class %in% is(object))
    invisible(object)
  else
    stop(paste('Object is not a member of class', class))

}

.GDALDataTypes <- c('Unknown', 'Byte', 'UInt16', 'Int16', 'UInt32',
                    'Int32', 'Float32', 'Float64', 'CInt16', 'CInt32',
                    'CFloat32', 'CFloat64')

setClass('GDALMajorObject',
         representation(handle = 'externalptr'))
 
getDescription <- function(object) {

  assertClass(object, 'GDALMajorObject')

  .Call('RGDAL_GetDescription', object, PACKAGE="rgdal")

}

getGDALVersionInfo <- function(str="--version") {
    stopifnot(is.character(str))
    .Call("RGDAL_GDALVersionInfo", str, PACKAGE="rgdal")
}

getGDALCheckVersion <- function() {
    .Call("RGDAL_GDALCheckVersion", PACKAGE="rgdal")
}

getGDAL_DATA_Path <- function() {
    res <- .Call("RGDAL_GDAL_DATA_Info", PACKAGE="rgdal")
    res <- sub("prime_meridian.csv", "", res)
    n <- nchar(res)
    res <- substring(res, 1, n-1)
    res
}

version_sp_linkingTo <- function() {
    .Call("rgdal_sp_linkingTo_version")
}




setClass('GDALDriver', 'GDALMajorObject')

setClass('GDALReadOnlyDataset', 'GDALMajorObject')

setClass('GDALDataset', 'GDALReadOnlyDataset')

setClass('GDALTransientDataset', 'GDALDataset')
         
setClass('GDALRasterBand', 'GDALMajorObject')

getGDALDriverNames <- function() {
  res <- .Call('RGDAL_GetDriverNames', PACKAGE="rgdal")
  has_isRaster <- 0L
  if (!is.null(attr(res, "isRaster"))) has_isRaster <- 1L
  if (has_isRaster) res$isRaster <- attr(res, "isRaster")
  res <- as.data.frame(res, stringsAsFactors=FALSE)
  if (has_isRaster) res <- res[res$isRaster,]
  res <- res[order(res$name),]
  row.names(res) <- NULL
  res
}

setMethod('initialize', 'GDALDriver',
          def = function(.Object, name, handle = NULL) {
            if (is.null(handle)) {
              slot(.Object, 'handle') <- {
                .Call('RGDAL_GetDriver', as.character(name), PACKAGE="rgdal")
              }
            } else {
              slot(.Object, 'handle') <- handle
            }
            .Object
          })

getDriverName <- function(driver) {

  assertClass(driver, 'GDALDriver')

  .Call('RGDAL_GetDriverShortName', driver, PACKAGE="rgdal")

}

getDriverLongName <- function(driver) {

  assertClass(driver, 'GDALDriver')

  .Call('RGDAL_GetDriverLongName', driver, PACKAGE="rgdal")

}

setMethod('initialize', 'GDALReadOnlyDataset',
          def = function(.Object, filename, silent=FALSE, handle = NULL) {
            if (is.null(handle)) {
              filename <- as.character(filename)
	      if (nchar(filename) == 0) stop("empty file name")
              silent <- as.logical(silent)
              if (length(silent) != 1L || is.na(silent) || !is.logical(silent))
                  stop("options(warn) not set")
              slot(.Object, 'handle') <- {
                .Call('RGDAL_OpenDataset',
                        normalizePath(filename, mustWork=FALSE), 
			TRUE, silent, PACKAGE="rgdal")
              }
            } else {
              slot(.Object, 'handle') <- handle
            }
            cfn <- function(handle) .Call('RGDAL_CloseHandle', 
		handle, PACKAGE="rgdal")
            .setCollectorFun(slot(.Object, 'handle'), cfn)
            .Object
          })

setMethod('initialize', 'GDALDataset',
          def = function(.Object, filename, silent=FALSE, handle = NULL) {
            if (is.null(handle)) {
              filename <- as.character(filename)
	      if (nchar(filename) == 0) stop("empty file name")
              silent <- as.logical(silent)
              if (length(silent) != 1L || is.na(silent) || !is.logical(silent))
                  stop("options(warn) not set")
              slot(.Object, 'handle') <- {
                .Call('RGDAL_OpenDataset', 
                        normalizePath(filename, mustWork=FALSE), 
			FALSE, silent, PACKAGE="rgdal")
              }
            } else {
              slot(.Object, 'handle') <- handle
            }
            cfn <- function(handle) .Call('RGDAL_CloseHandle', 
		handle, PACKAGE="rgdal")
            .setCollectorFun(slot(.Object, 'handle'), cfn)
            .Object
          })

setMethod('initialize', 'GDALTransientDataset',
          def = function(.Object, driver, rows, cols, bands = 1,
            type = 'Byte', options = NULL, fname = NULL, handle = NULL) {
            if (is.null(handle)) {
              typeNum <- match(type, .GDALDataTypes, 1) - 1
	      if (is.null(fname)) {
                  my_tempfile <- tempfile()
              } else {
                  my_tempfile <- paste(tempdir(), "/",
                      paste(sample(letters, 3), collapse=""),
                      basename(fname[1]), sep="")
              }
	      if (nchar(my_tempfile) == 0) stop("empty file name")
	      if (!is.null(options)) options <- as.character(options)
              slot(.Object, 'handle') <- .Call('RGDAL_CreateDataset', driver,
                                              as.integer(c(cols, rows, bands)),
                                              as.integer(typeNum),
                                              options,
                                              my_tempfile, PACKAGE="rgdal")
            } else {
              slot(.Object, 'handle') <- handle
            }
            cfn <- function(handle) .Call('RGDAL_CloseHandle',
#            cfn <- function(handle) .Call('RGDAL_CloseDataset', RSB 081030
		handle, PACKAGE="rgdal")
            .setCollectorFun(slot(.Object, 'handle'), cfn)
            .Object
          })

getDriver <- function(dataset) {

  assertClass(dataset, 'GDALReadOnlyDataset')

  new('GDALDriver',
      handle = .Call('RGDAL_GetDatasetDriver', dataset, PACKAGE="rgdal"))

}

copyDataset <- function(dataset, driver, strict = FALSE, options = NULL, fname = NULL) {

  assertClass(dataset, 'GDALReadOnlyDataset')
  
  if (missing(driver)) driver <- getDriver(dataset)
  else if (is.character(driver)) driver <- new("GDALDriver", driver)

  if (is.null(fname)) {
     my_tempfile <- tempfile()
  } else {
     my_tempfile <- paste(tempdir(), "/",
       paste(sample(letters, 3), collapse=""),
       basename(fname[1]), sep="")
  }
  if (nchar(my_tempfile) == 0) stop("empty file name")
  if (!is.null(options) && !is.character(options))
    stop("options not character")
  
  new.obj <- new('GDALTransientDataset',
                 handle = .Call('RGDAL_CopyDataset',
                   dataset, driver,
                   as.integer(strict),
                   as.character(options),
                   my_tempfile, PACKAGE="rgdal"))

  new.obj
  
}

saveDataset <- function(dataset, filename, options=NULL, returnNewObj=FALSE) {

  assertClass(dataset, 'GDALReadOnlyDataset')
  
  new.class <- ifelse(class(dataset) == 'GDALTransientDataset',
                      'GDALDataset', class(dataset))
  if (!is.null(options) && !is.character(options))
    stop("options not character")
  
  filename <- as.character(filename)
  if (nchar(filename) == 0) stop("empty file name")
  new.obj <- new(new.class,
                 handle = .Call('RGDAL_CopyDataset',
                   dataset, getDriver(dataset),
                   FALSE, options, normalizePath(filename, mustWork=FALSE),
                   PACKAGE="rgdal"))

  if (returnNewObj) return(new.obj)
  invisible(GDAL.close(new.obj))
}

setGeneric('closeDataset', function(dataset) standardGeneric('closeDataset'))

"closeDataset.default" <- function(dataset) 
	stop("No default method for closeDataset")

setMethod("closeDataset", signature("ANY"), closeDataset.default)

setMethod('closeDataset', 'GDALReadOnlyDataset',
          def = function(dataset) {
            .setCollectorFun(slot(dataset, 'handle'), NULL)
            .Call('RGDAL_CloseDataset', dataset, PACKAGE="rgdal")
            invisible(gc())
          })

setMethod('closeDataset', 'GDALTransientDataset',
          def = function(dataset) {
            driver <- getDriver(dataset)
#            filename <- getDescription(dataset)
#            .Call('RGDAL_CloseDataset', driver, filename, PACKAGE="rgdal")
            .Call('RGDAL_CloseDataset', driver, PACKAGE="rgdal")
            invisible(gc())
            callNextMethod()
          })


saveDatasetAs <- function(dataset, filename, driver = NULL, options=NULL) {

  .Deprecated("saveDataset")

  assertClass(dataset, 'GDALReadOnlyDataset')
  
  filename <- as.character(filename)
  if (nchar(filename) == 0) stop("empty file name")
  if (is.null(driver)) driver <- getDriver(dataset)
  if (!is.null(options) && !is.character(options))
    stop("options not character")
  
  new.obj <- new('GDALReadOnlyDataset',
                 handle = .Call('RGDAL_CopyDataset',
                   dataset, driver, FALSE, options,
                   normalizePath(filename, mustWork=FALSE), PACKAGE="rgdal"))
  
  closeDataset(new.obj)
  
  err.opt <- getOption('show.error.messages')

  options(show.error.messages = FALSE)

  new.obj <- try(new('GDALDataset', filename))

  options(show.error.messages = err.opt)

  if (inherits(new.obj, 'try-error'))
    new.obj <- new('GDALReadOnlyDataset', filename)

  closeDataset(dataset)

  eval.parent(dataset <- new.obj)

  invisible(new.obj)
  
}


deleteDataset <- function(dataset) {

  assertClass(dataset, 'GDALDataset')
  
  driver <- getDriver(dataset)
  
  filename <- getDescription(dataset)
  
  .Call('RGDAL_DeleteFile', driver, filename, PACKAGE="rgdal")
  
  closeDataset(dataset)

}

isObjPtrNULL <- function(ptr) {

  stopifnot(is(ptr, "GDALMajorObject"))

  .Call("isGDALObjPtrNULL", ptr, PACKAGE="rgdal")

}

GDAL.open <- function(filename, read.only = TRUE, silent = FALSE) {
  
	res <- if(read.only)
          new("GDALReadOnlyDataset", filename, silent=silent)
        else
          new("GDALDataset", filename, silent=silent)
        
	res
        
}

GDAL.close <- function(dataset) {
            isTrans <- is(dataset, "GDALTransientDataset")
            if (isTrans) {
                if (isObjPtrNULL(dataset)) stop("object already closed")
#                filename <- getDescription(dataset)
            }
            .setCollectorFun(slot(dataset, 'handle'), NULL)
            .Call('RGDAL_CloseDataset', dataset, PACKAGE="rgdal")
#            .Call("RGDAL_CloseHandle", slot(dataset, 'handle'),
#                PACKAGE="rgdal")
             invisible(NULL)
#            if (isTrans) {
#                basen <- basename(filename)
#                dirn <- dirname(filename)
#                lf <- list.files(path=dirn, pattern=basen)
#                flf <- paste(dirn, lf, sep="/")
#                unlink(flf)
#            }
#            invisible(gc())
}

setMethod('dim', 'GDALReadOnlyDataset',
          def = function(x) {
            nrows <- .Call('RGDAL_GetRasterYSize', x, PACKAGE="rgdal")
            ncols <- .Call('RGDAL_GetRasterXSize', x, PACKAGE="rgdal")
            nbands <- .Call('RGDAL_GetRasterCount', x, PACKAGE="rgdal")
            if (nbands < 1) warning("no bands in dataset")
            if (nbands > 1)
              c(nrows, ncols, nbands)
            else
              c(nrows, ncols)
          })

getProjectionRef <- function(dataset, OVERRIDE_PROJ_DATUM_WITH_TOWGS84=NULL) {

  assertClass(dataset, 'GDALReadOnlyDataset')

  vs <- strsplit(strsplit(getGDALVersionInfo(), ",")[[1]][1], " ")[[1]][2]
  env_absent <- is.null(getCPLConfigOption("OVERRIDE_PROJ_DATUM_WITH_TOWGS84"))
  if ((vs > "1.8.0") && env_absent) {
    if (is.null(OVERRIDE_PROJ_DATUM_WITH_TOWGS84))
      OVERRIDE_PROJ_DATUM_WITH_TOWGS84 <- get_OVERRIDE_PROJ_DATUM_WITH_TOWGS84()
    stopifnot(is.logical(OVERRIDE_PROJ_DATUM_WITH_TOWGS84))
    stopifnot(length(OVERRIDE_PROJ_DATUM_WITH_TOWGS84) == 1)
    if (!OVERRIDE_PROJ_DATUM_WITH_TOWGS84) {
      setCPLConfigOption("OVERRIDE_PROJ_DATUM_WITH_TOWGS84", "NO")
      res <- .Call('RGDAL_GetProjectionRef', dataset, PACKAGE="rgdal")
      setCPLConfigOption("OVERRIDE_PROJ_DATUM_WITH_TOWGS84", NULL)
    } else {
      res <- .Call('RGDAL_GetProjectionRef', dataset, PACKAGE="rgdal")
    }
  } else {
    res <- .Call('RGDAL_GetProjectionRef', dataset, PACKAGE="rgdal")
  }
  res
}

putRasterData <- function(dataset,
                          rasterData,
                          band = 1,
                          offset = c(0, 0)) {

  assertClass(dataset, 'GDALDataset')

  offset <- rep(offset, length.out = 2)
  
  raster <- getRasterBand(dataset, band)
  
  .Call('RGDAL_PutRasterData', raster, rasterData, 
	as.integer(offset), PACKAGE="rgdal")

}

getRasterTable <- function(dataset,
                           band = NULL,
                           offset = c(0, 0),
                           region.dim = dim(dataset)) {

  assertClass(dataset, 'GDALReadOnlyDataset')

  offset <- rep(offset, length.out = 2)
  region.dim <- rep(region.dim, length.out = 2)

  rasterData <- getRasterData(dataset, band,
                              offset = offset,
                              region.dim = region.dim, list_out=TRUE)

  if (is.null(band)) {

    nbands <- .Call('RGDAL_GetRasterCount', dataset, PACKAGE="rgdal")
    if (nbands < 1) stop("no bands in dataset")
    band <- 1:nbands

  } else {

    nbands <- length(band)

  }

#  dim(rasterData) <- c(region.dim, nbands)

  geoTrans <- getGeoTransFunc(dataset)

  y.i <- 1:region.dim[1] - 0.5 + offset[1]
  x.i <- 1:region.dim[2] - 0.5 + offset[2]

  y.i <- rep(y.i, each = length(x.i))
  x.i <- rep(x.i, len = prod(region.dim))

  out <- geoTrans(x.i, y.i)

#  out <- cbind(out$x, out$y)
  out <- data.frame(x=out$x, y=out$y)
  rasterData <- as.data.frame(rasterData)

#  for (b in band) { 
#    vec <- as.numeric(rasterData[, , b])
#    out <- cbind(out, vec)
#  }

#  out <- as.data.frame(out)
    
#  names(out) <- c('x', 'y', paste('band', 1:nbands, sep = ''))
  out <- cbind(out, rasterData)

  out

}
                           
getRasterData <- function(dataset,
                          band = NULL,
                          offset = c(0, 0),
                          region.dim = dim(dataset),
                          output.dim = region.dim,
                          interleave = c(0, 0),
                          as.is = FALSE, list_out=FALSE) {

    assertClass(dataset, 'GDALReadOnlyDataset')

    offset <- rep(offset, length.out = 2)
    region.dim <- rep(region.dim, length.out = 2)
    output.dim <- rep(output.dim, length.out = 2)
    interleave <- rep(interleave, length.out = 2)

    nbands <- .Call('RGDAL_GetRasterCount', dataset, PACKAGE="rgdal")
    if (nbands < 1) stop("no bands in dataset")

    if (is.null(band)) band <- 1:nbands
  
    x <- array(dim = as.integer(c(rev(output.dim), length(band))))
    for (i in seq(along = band)) {

        raster <- getRasterBand(dataset, band[i])

        x[,,i] <- .Call('RGDAL_GetRasterData', raster,
                      as.integer(c(offset, region.dim)),
                      as.integer(output.dim),
                      as.integer(interleave),
                      PACKAGE="rgdal")
  
    }
    if (!as.is) {
        for (i in seq(along = band)) {

            raster <- getRasterBand(dataset, band[i])
            scale <- .Call('RGDAL_GetScale', raster, PACKAGE="rgdal")
            offset <- .Call('RGDAL_GetOffset', raster, PACKAGE="rgdal")

            if (scale != 1) x[,,i] <- x[,,i] * scale
            if (offset != 0) x[,,i] <- x[,,i] + offset
        }
    }
    if (!list_out) {
        if (length(band) == 1L) x <- drop(x)
        return(x)
    } else {
        X <- vector(mode="list", length=length(band))
        names(X) <- paste("band", 1:length(band), sep="")

        for (i in seq(along = band)) {

            X[[i]] <- as.vector(x[,,i])

            if (!as.is) {
  
                raster <- getRasterBand(dataset, band[i])
    
                catNames <- .Call('RGDAL_GetCategoryNames', raster,
                    PACKAGE="rgdal")
  
                if (!is.null(catNames)) {
                    ux <- sort(unique(na.omit(X[[i]])))
                    lCN <- length(catNames)
                    levels <- ((1:lCN)-1)
                    back_incls <- ux %in% levels
                    if (all(back_incls)) {
                        X[[i]] <- factor(X[[i]], levels=levels, labels=catNames)
                        if (!get("silent", envir=.RGDAL_CACHE)) {
                            cat("Input level values and names\n")
                            cat(paste(levels, " ", catNames, "\n", sep=""),
                                sep="")
                        }
                    } else {
                        warning("Assign CategoryNames manually, level/label mismatch")
                    }
                }
            }
        }
        return(X)
    }
}

getCategoryNames <- function(dataset, band = 1) {

  assertClass(dataset, 'GDALReadOnlyDataset')

  raster <- getRasterBand(dataset, band)
  
  catNames <- .Call('RGDAL_GetCategoryNames', raster, PACKAGE="rgdal")

  catNames
}

getBandColorTable <- function(raster) {

  ctab <- .Call('RGDAL_GetColorTable', raster, PACKAGE="rgdal") / 255

  if (length(ctab) == 0L) return(NULL)

  if (.Call('RGDAL_GetColorInterp', raster, PACKAGE="rgdal") == 'Palette')
    switch(.Call('RGDAL_GetPaletteInterp', raster, PACKAGE="rgdal"),  
           RGB = rgb(ctab[,1], ctab[,2], ctab[,3]),
           HSV = hsv(ctab[,1], ctab[,2], ctab[,3]), # Doesn't actually exist
           Gray = gray(ctab[,1]),
           gray(apply(ctab, 2, mean)))
  else
    gray(ctab[,1])

}

getColorTable <- function(dataset, band = 1) {

  assertClass(dataset, 'GDALReadOnlyDataset')

  raster <- getRasterBand(dataset, band)
  
  getBandColorTable(raster)
}

RGB2PCT <- function(x, band, driver.name = 'MEM',
                    ncolors = 256, set.ctab = TRUE) {
  
  assertClass(x, 'GDALReadOnlyDataset')

  if (ncolors < 2 || ncolors > 256)
    stop('Number of colors must be between 2 and 256')

  band <- rep(band, length.out = 3)

  dithered <- new('GDALTransientDataset',
                  new('GDALDriver', driver.name),
                  nrow(x), ncol(x))

  ctab <- .Call('RGDAL_GenCMap',
                getRasterBand(x, band[1]),
                getRasterBand(x, band[2]),
                getRasterBand(x, band[3]),
                getRasterBand(dithered),
                as.integer(ncolors),
                as.logical(set.ctab),
                PACKAGE = "rgdal") / 255

  if (set.ctab)
    dithered
  else
    list(dataset = dithered,
         pct = rgb(ctab[,1], ctab[,2], ctab[,3]))

}  

displayDataset <- function(x, offset = c(0, 0), region.dim = dim(x),
                           reduction = 1, band = 1, col = NULL,
                           reset.par = TRUE, max.dim = 500, ...) {

  assertClass(x, 'GDALReadOnlyDataset')

  offset <- rep(offset, length.out = 2)
  region.dim <- rep(region.dim, length.out = 2)
  reduction <- rep(reduction, length.out = 2)

  offset <- offset %% dim(x)[1:2]
  
  oob <- (region.dim + offset) > dim(x)[1:2]
  
  if (any(oob)) region.dim[oob]  <-  dim(x)[oob] - offset[oob]

  reduction[reduction < 1] <- 1

  plot.dim <- region.dim / reduction
            
  if (any(plot.dim > max.dim))
    plot.dim <- max.dim * plot.dim / max(plot.dim)

  image.data <- getRasterData(x, band[1], offset,
                              region.dim, plot.dim,
                              as.is = TRUE)
#  image.data <- array(image.data[[1]], t(plot.dim))

  if (is.null(col)) {
    
    max.val <- max(image.data, na.rm = TRUE)

    if (!is.finite(max.val)) {
      image.data[] <- 2
      max.val <- 2
    }

    col <- getColorTable(x, band)[1:(max.val + 1)]
      
  }

  if (is.null(col)) col <- gray(seq(0, 1, len = 256))
  
  par.in <- par(no.readonly = TRUE)

  if (reset.par) on.exit(par(par.in))

  par(pin = max(par.in$pin)
      * par.in$fin / max(par.in$fin)
      * rev(plot.dim) / max(plot.dim))
  
  image.data <- image.data[, ncol(image.data):1]

  image.default(image.data + 1, col = col, axes = FALSE, ...)
            
  invisible(list(image.data = image.data, col = col, par.in = par.in))

}

setMethod('initialize', 'GDALRasterBand',
          def =  function(.Object, dataset, band = 1) {
            slot(.Object, 'handle') <- .Call('RGDAL_GetRasterBand',
                                            dataset, as.integer(band), 
					    PACKAGE="rgdal")
            .Object
          })

setMethod('dim', 'GDALRasterBand',
          def = function(x) {
            c(.Call('RGDAL_GetYSize', x, PACKAGE="rgdal"),
              .Call('RGDAL_GetXSize', x, PACKAGE="rgdal"))
          })

getGeoTransFunc <- function(dataset) {

  assertClass(dataset, 'GDALReadOnlyDataset')

  geoTrans <- .Call('RGDAL_GetGeoTransform', dataset, PACKAGE="rgdal")
  if (attr(geoTrans, "CE_Failure")) warning("GeoTransform values not available")
  rotMat <- matrix(geoTrans[c(2, 3, 5, 6)], 2)

  offset <- geoTrans[c(1, 4)]

  function(x, y) {

    x <- cbind(x, y)

    x <- x %*% rotMat

    list(x = x[,1] + offset[1],
         y = x[,2] + offset[2])

  }

}

getRasterBand <- function(dataset, band = 1) {

  assertClass(dataset, 'GDALReadOnlyDataset')

  new('GDALRasterBand', dataset, band)

}

getRasterBlockSize <- function(raster) {

  assertClass(raster, 'GDALRasterBand')
  
  .Call('RGDAL_GetRasterBlockSize', raster, PACKAGE="rgdal")
  
}

get_OVERRIDE_PROJ_DATUM_WITH_TOWGS84 <- function() {
  get("OVERRIDE_PROJ_DATUM_WITH_TOWGS84", envir=.RGDAL_CACHE)
}

set_OVERRIDE_PROJ_DATUM_WITH_TOWGS84 <- function(value) {
        stopifnot(is.logical(value))
        stopifnot(length(value) == 1)
        assign("OVERRIDE_PROJ_DATUM_WITH_TOWGS84", value, envir = .RGDAL_CACHE)
        get_OVERRIDE_PROJ_DATUM_WITH_TOWGS84()
}

getCPLConfigOption <- function(ConfigOption) {
    stopifnot(is.character(ConfigOption))
    stopifnot(length(ConfigOption) == 1)
    .Call("RGDAL_CPLGetConfigOption", ConfigOption, PACKAGE="rgdal")
}

setCPLConfigOption <- function(ConfigOption, value) {
    stopifnot(is.character(ConfigOption))
    stopifnot(length(ConfigOption) == 1)
    if (!is.null(value)) {
        stopifnot(is.character(value))
        stopifnot(length(value) == 1)
    }
    .Call("RGDAL_CPLSetConfigOption", ConfigOption, value, PACKAGE="rgdal")
    .Call("RGDAL_CPLGetConfigOption", ConfigOption, PACKAGE="rgdal")
}

GDAL_iconv <- function() {
    .Call("RGDAL_CPL_RECODE_ICONV", PACKAGE="rgdal")
}

