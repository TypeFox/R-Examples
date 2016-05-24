get.sst <-
function(xlim = NULL, ylim = NULL, tlim = NULL, server = "http://atlas.nmfs.hawaii.edu/cgi-bin/reynolds_extract.py") {
    ## require(sp)  ## removed with new namespace import MDS 2011-10-06
    ## TODO:
    ## checks on coordinate ranges, lengths
    ## file handling, to cache them somehow
    ## checks on file sizes
    if (is.null(xlim)) xlim <- c(-180, 180)
    if (is.null(ylim)) ylim  <- c(-90, 90)
    ## work in POSIXct
    if(is.null(tlim)) {
        tlim <- rep(Sys.time(), 2)
    }

    tlim <- tlim + c(-1, 0) * 5* 24 * 3600
    yrRan <- as.POSIXlt(tlim)$year + 1900
    daysIntoYearRan <- as.POSIXlt(tlim)$yday + 1
    tfile <- tempfile()
    string <- ""
    string <- paste(string, "?lon1=", xlim[1], "&lon2=", xlim[2],
                    string, "&lat1=", ylim[1], "&lat2=", ylim[2],
                    string, "&year1=", yrRan[1], "&day1=", daysIntoYearRan[1],
                    string, "&year2=", yrRan[2], "&day2=", daysIntoYearRan[2],
                    sep = "")
    link <- file.path(server, string)
    download.file(link, tfile, mode = "wb")
    files <- unzip(tfile, NULL, exdir =dirname(tfile))

    files0 <- files
    ## need more cases than this (expected minimal file size? files equal size within some tolerance?)
    ## try/catches for conversion to gridded?
    files <- files[file.info(files)$size > 0]

    times <-  as.POSIXct(strptime(basename(files), "RS%Y%j"))
    ord <- order(times)
    files <- files[ord]
    times <- times[ord]


    if (length(files) < 1) stop("no valid files downloaded", files)
    grid0 <- read.table(files[1])

    coordinates(grid0) <- 2:1
    gridded(grid0) <- TRUE
    grid0 <- as.image.SpatialGridDataFrame(grid0)
    grid0$x <- grid0$x - 360
    grid0$t <- times

   if (length(files) > 1) {
        grid0$z <- array(grid0$z, c(dim(grid0$z), length(files)))
        for (i in 2:length(files)) {
          gridi <- read.table(files[i])

          coordinates(gridi) <- 2:1
          gridded(gridi) <- TRUE
          gridi <- as.image.SpatialGridDataFrame(gridi)
          grid0$z[,,i] <- gridi$z

      }
    }

    #class(grid0) <- c("genimg", class(grid0))
    grid0
}

#get.reyn.sst <- function(filename = "G:/DATA/Reynolds/sst.wkmean.1990-present.nc") {
#    stopifnot(file.exists(filename))
#    stopifnot(require(RNetCDF))

#    con <- open.nc(filename)
#    lon <-


#mds}
