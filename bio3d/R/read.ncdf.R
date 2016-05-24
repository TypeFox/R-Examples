`read.ncdf` <-
function(trjfile, headonly = FALSE, verbose = TRUE, time=FALSE,
         first = NULL, last= NULL, stride = 1, cell = FALSE,
         at.sel = NULL){
  # Adding option 'at.sel' to select a subset of the structure using
  # an object of class 'select'
 
  # Currently file open in SHARE mode is not supported by NCDF package.
  # Multicore support is suppressed
#  ncore = 1

  ##- Load ncdf4 package
  oops <- requireNamespace("ncdf4", quietly = TRUE)
  if(!oops)
    stop("Please install the ncdf4 package from CRAN")
  
  # open files 
  nc <- lapply(trjfile, function (fnm) {

     nc <- ncdf4::nc_open(fnm, readunlim=FALSE)
     conv <- ncdf4::ncatt_get( nc, varid=0, "Conventions")$value
     if(conv!="AMBER")
        warning(paste("File conventions is not set to AMBER",
                  fnm, sep=": "))
     return(nc)
  })

  # set first and last frame no for each file
  # used for time/frame range selection
  flen <- unlist(lapply(nc, function(n) return(n$dim$frame$len)))
  frange <- matrix(c(1, cumsum(flen[-length(nc)])+1, cumsum(flen)), 
                 nrow=length(nc))

  if(verbose) {
    if(length(trjfile)>1)
       print(paste("Reading ", length(trjfile), "files"))
    else
       print(paste("Reading file", trjfile))
    print(paste("Produced by program:",
      ncdf4::ncatt_get( nc[[1]], varid=0, "program")$value))
    print(paste("File conventions",
      ncdf4::ncatt_get( nc[[1]], varid=0, "Conventions")$value, "version",
      ncdf4::ncatt_get( nc[[1]], varid=0, "ConventionVersion")$value))
    print(paste("Frames:", sum(flen)))
    print(paste("Atoms:", nc[[1]]$dim$atom$len))
  }
 
  # read heads, cell, or coordinates
  retval <- lapply(seq_along(nc), function (inc) {
    first.atom <-  1
    count.atom <- -1
    if(!is.null(at.sel)) {
      if(!is.select(at.sel)) stop("'at.sel' must be an object of class 'select'. See 'atom.select'.")
      atom.ind <- xyz2atom(at.sel$xyz)
      first.atom <- min(atom.ind)
      count.atom <- diff(range(atom.ind)) + 1
    }
    
     nc <- nc[[inc]]  
     frange <- frange[inc,]
     ss = 1
     ee = nc$dim$frame$len
     if (!is.null(c(first, last))) {

        #check frame No or time
        if(time) {
           btime <- ncdf4::ncvar_get(nc, "time", 1, 1)
           etime <- ncdf4::ncvar_get(nc, "time", nc$dim$frame$len, 1)
        } else {
           btime = frange[1]
           etime = frange[2]
        }

        if((!is.null(first) && (etime < first)) || 
              (!is.null(last) && last >=0 && btime > last) ||
             (!is.null(first) && !is.null(last) && last >=0 && last < first) ) {
#           if(verbose) print(paste("Skip file", nc$filename))
           ncdf4::nc_close(nc)
           return()
        } 
        if(!headonly) {
           timeall <- btime:etime
           if(time) timeall <- ncdf4::ncvar_get(nc, "time")

           ss <- if(is.null(first)) 1 else which((timeall - first) >=0 )[1]
           if(is.null(last) || last < 0 || last > etime) {
              ee = nc$dim$frame$len
           } else {
              ee <- which((timeall - last) <= 0) 
              ee <- ee[length(ee)]
           }
        }
     }
     tlen = ee - ss + 1
     conv <- ncdf4::ncatt_get( nc, varid=0, "Conventions")$value
     if(headonly) {
       ## Report only header information
       return(list("file"=nc$filename,
                   "conv"=conv,
                   "frames"=nc$dim$frame$len,
                   "atoms"=nc$dim$atom$len))
     ##time  <- ncdf4::ncvar_get(nc,"time")
     }
     if(cell) {
        celldata <- ncdf4::ncvar_get(nc, "cell_lengths", c(1, ss), 
                             c(-1, tlen))
        celldata <- t( rbind(celldata, ncdf4::ncvar_get(nc, "cell_angles", 
                         c(1, ss), c(-1, tlen))) )
        if(time)
           rownames(celldata) <- ncdf4::ncvar_get(nc, "time", ss, tlen)
        ncdf4::nc_close(nc)
        return( celldata )
     }
     if(count.atom < 0) count.atom = nc$dim$atom$len

     # To solve 32-bit limitation problem with large trajectory file
     .check <- (3 * count.atom * tlen) / (2^31 - 1)
     if(.check > 1) {
        .nb <- floor(.check) + 1
        .nn <- floor(tlen / .nb)
        .ss <- seq(ss, ss + tlen - 1, .nn)
        .tlen <- rep(.nn, length(.ss) - 1)
        .tlen <- c(.tlen, tlen - sum(.tlen))
        coords <- sapply(1:length(.ss), function(i) 
            ncdf4::ncvar_get(nc, "coordinates", c(1, first.atom, .ss[i]), 
                          c(-1, count.atom, .tlen[i])) )
     } else {    
        coords <- ncdf4::ncvar_get(nc, "coordinates", c(1, first.atom, ss), 
                          c(-1, count.atom, tlen))
     }
     if(!is.null(at.sel)) coords <- coords[,atom.ind - first.atom + 1,]
     coords <- matrix( coords, ncol=(dim(coords)[2]*3), byrow=TRUE )
     if(time)
        rownames(coords) <- ncdf4::ncvar_get(nc, "time", ss, tlen)
     ncdf4::nc_close(nc)
     return( coords )
  } )

  if(headonly) {
    retval <- do.call("c", retval)
  }
  else if(cell) {
    retval <- do.call(rbind, retval)
    ##retval <- as.data.frame(retval, stringsAsFactors=FALSE)
  }
  else {
    retval <- do.call(rbind, retval)

    ## take every "stride" frame
    retval <- as.xyz(subset(retval, (1:nrow(retval)) %in% seq(1, nrow(retval), stride)))
  }

  return( retval )
}
