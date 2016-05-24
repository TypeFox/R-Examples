
#' Determine or modify file name extensions
#' 
#' Function \code{get.file.extension} determines the file extension, \code{set.file.extension} changes it, and \code{default.file.extension} changes it only if it is not already specified.
#' @name set.file.extension
#' @param filename character vector: file name(s), possibly including paths and extensions; a file name ending with a \code{"."} is interpreted as having extension \code{""}, while a file name that doesn't contain a \code{"."} is interpreted has having no extension.
#' @param extension character string: file extension, without the dot
#' @param fsep character: separator between paths
#' @param force logical argument to \code{default.file.extension}: force the file extension to be \code{extension} (same result as \code{set.file.extension}), or only set it to extension if it has not been specified?
#' @return character vector of same length as \code{filename}
#' @examples 
#' fnm = c("C:/TEMP.DIR/temp","C:/TEMP.DIR/tmp.txt","tempfile.")
#' get.file.extension(fnm)
#' set.file.extension(fnm,extension=".TMP")
#' default.file.extension(fnm,extension=".TMP")
#' @keywords file utilities
#' @export
set.file.extension = function(filename, extension, fsep=.Platform$file.sep) {
    if (extension=="") extension = "."
    if (substr(extension,1,1)!=".") extension = paste(".",extension,sep="")
    if (Sys.info()["sysname"] == "Windows")
        filename = gsub("\\",fsep,filename,fixed=TRUE)
    the.extension = get.file.extension(filename)
    filename = substr(filename,1,nchar(filename)-nchar(the.extension))
    return( paste( filename, extension, sep="") )
}


#' @rdname set.file.extension
#' @name get.file.extension
#' @export
get.file.extension = function(filename, fsep=.Platform$file.sep) {
    ext = rep("",length(filename))
    has.dot.extension = substring(filename, nchar(filename))=="."
    if (Sys.info()["sysname"] == "Windows")
        filename = gsub("\\",fsep,filename,fixed=TRUE)
    split = strsplit(filename,fsep,fixed=TRUE)
    split = sapply( split, function(x) x[length(x)] )
    split = strsplit(split,".",fixed=TRUE)
    ext = sapply( split, function(x) x[length(x)] )
    has.extension = sapply(split,length) > 1
    ext[ has.extension ] = paste(".",ext[has.extension],sep="")
    ext[ !has.extension ] = ""
    ext[ has.dot.extension ] = "."
    return(ext)
}


#' @rdname set.file.extension
#' @name default.file.extension
#' @export
default.file.extension = function(filename, extension, force=FALSE) {
    if (force) {
        filename = set.file.extension(filename,extension)
    } else {
        use.default = (get.file.extension(filename)=="")
        if (any(use.default))
            filename[use.default] = set.file.extension(filename[use.default],extension)
    }
    return(filename)
}


#' Convert file name to variable name
#'
#' Convert a file name into a variable name
#' @name create.variable.name
#' @param filename character string
#' @param prefix character string: optional prefix to be added
#' @param fsep character used to separate path components
#' @examples 
#' \dontrun{
#' create.variable.name("C:/my-path/my-file-name.Rd",prefix="res")
#' }
#' @keywords utilities
#' @export
create.variable.name = function( filename, prefix = NULL, fsep = .Platform$file.sep )
{
    has.dot.extension = substring(filename, nchar(filename))=="."
    varname = filename
    if (Sys.info()["sysname"] == "Windows")
        varname = gsub("\\",fsep,varname,fixed=TRUE)
    varname = strsplit(varname,fsep,fixed=TRUE)[[1]]
    varname = varname[length(varname)]
    varname = strsplit(varname,".",fixed=TRUE)[[1]]
    if (length(varname) > 1) {
        if (!has.dot.extension)
            varname = varname[ 1 : (length(varname)-1) ]
        varname = paste(varname,collapse=".")
    }
    varname = gsub("-",".",gsub("_",".",varname))
    if (!is.null(prefix)) if (prefix!="") 
        varname = paste(prefix,".",varname,sep="")
    return(varname)
}



#' Read/write ASCII, SAGA and Rd Grid Files
#'
#' These functions provide simple interfaces for reading and writing grids from/to ASCII grids and Rd files. Grids are stored as matrices, their headers in lists.
#' @name read.ascii.grid
#' @param file file name of an ASCII grid (extension defaults to \code{.asc} if not specified), or a connection open for reading or writing, as required
#' @param fname file name of a grid stored as an R (\code{.Rd}) file; extension defaults to \code{.Rd}
#' @param return.header logical: should the grid header be returned (default), or just the grid data matrix? In the former case, \code{read.ascii.grid} returns a list with two components named \code{data} and \code{header}.
#' @param print numeric, specifying how detailed the output reporting the progress should be (currently 0 to 2, 0 being minimum output).
#' @param nodata.values optional numeric vector specifying nodata values to be used in addition to the nodata value specified in the grid header; nodata values are converted to \code{NA}.
#' @param at.once logical: if \code{TRUE}, read the whole grid with one \code{scan} command; if \code{FALSE}, read it row by row using \code{scan} with option \code{nlines=1}.
#' @param data grid data: a data matrix, or a list with components \code{data} (the grid data matrix) and \code{header} (the grid header information).
#' @param header optional list argument specifying the grid header information as returned by the \code{read.ascii.grid} or \code{read.ascii.grid.header} function; see Details
#' @param write.header logical: should the header be written with the grid data? (default: \code{TRUE}) 
#' @param digits numeric: if not missing, write data rounded to this many decimal places
#' @param hdr.digits numeric: see \code{hdr.prec}
#' @param hdr.prec numeric: write (non-integer) header data with this many decimal places; a value of 9 or higher is recommended for compatibility with SAGA GIS (default: 10)
#' @param dec character (default: \code{"."}): decimal mark used in input or output file
#' @param georef character: specifies whether the output grid should be  georeferenced by the \code{"center"} or \code{"corner"} of its lower left grid cell; defaults to \code{"corner"}.
#' @param compress logical: should the \code{.Rd} file written by \code{write.Rd.file} be compressed? (default: \code{TRUE})
#' @param prec integer: number of digits of temporary ASCII grid used for importing or exporting a SAGA grid
#' @param na.strings passed on to \code{\link{scan}}.
#' @param ... \code{read.sgrd}, \code{write.sgrd}: additional arguments to be passed to \code{rsaga.geoprocessor}
#'
#' @return The \code{read.*} functions return either a list with components \code{data} (the grid data matrix) and \code{header} (the grid header  information, see below), if \code{return.header=TRUE}, or otherwise  just the grid data matrix \code{return.header=FALSE}.
#' 
#' The grid data matrix is a numeric matrix whose first column corrensponds to the first (i.e. northernmost) row of the grid. Columns run from left = West to right = East.
#' 
#' The header information returned by the \code{read.ascii.grid[.header]} functions (if \code{return.header=TRUE}) is a list with the following components:
#'  \item{ncols}{Number of grid columns.}
#'  \item{nrows}{Number of grid rows.}
#'  \item{xllcorner}{x coordinate of the corner of the lower left grid cell.}
#'  \item{yllcorner}{y coordinate of the corner of the lower left grid cell.}
#'  \item{cellsize}{Single numeric value specifying the size of a grid cell or pixel in both x and y direction.}
#'  \item{nodata_value}{Single numeric value being interpreted as \code{NA} (typically \code{-9999}.}
#'  \item{xllcenter}{x coordinate of the center of the lower left grid cell}
#'  \item{yllcenter}{y coordinate of the center of the lower left grid cell}
#' Note: The order of the components, especially of \code{?llcorner} and \code{?llcenter}, may change, depending on the order in which they appear in the grid header and on the georeferencing method (center or corner) used for the grid. The \code{?llcorner} and \code{?llcenter} attributes differ only by \code{cellsize/2}.
#' @author Alexander Brenning
#' @note \code{read.sgrd} and \code{write.sgrd} import/export grids indirectly by creating temporary ASCII grid files (this explains why \code{write.sgrd} has \code{prec} and \code{hdr.prec} arguments). Consider using \code{readGDAL} and \code{\link[rgdal]{writeGDAL}} in package \code{rgdal} instead, which are likely more efficient but may require coercion of your gridded data to/from a \code{Spatial...DataFrame-class}.
#'
#' The \code{read.Rd.grid} and \code{write.Rd.grid} functions use the \code{load} and \code{save} commands to store a grid. The variable name used is \code{data}, which is either a numeric matrix or a list with components \code{data} (the grid data matrix) and \code{header} (the grid header information).
#' @seealso \code{readGDAL} and \code{\link[rgdal]{writeGDAL}} in package \code{rgdal}, and \code{readAsciiGrid} and \code{\link[maptools]{writeAsciiGrid}} in package \code{maptools}
#' @keywords file spatial interface
#' @export
read.ascii.grid = function( file, return.header = TRUE, print = 0,
    nodata.values = c(), at.once = TRUE, na.strings = "NA" )
{
    if (is.character(file)) {
        file = default.file.extension(file,".asc")
        con = file(file,open="r")
        on.exit(close(con), add = TRUE)
    } else {
        con = file # was missing - bug fixed 2008-05-03
        if (!isOpen(file,"read"))
            stop("'file' must be a file name or a connection opened for reading")
    }
    hdr = read.ascii.grid.header(con)
    if (at.once) {
        data = scan(con, nlines=hdr$nrows, quiet=TRUE, na.strings = na.strings)
        data = matrix(data, ncol=hdr$ncols, nrow=hdr$nrows, byrow=TRUE)
        na = is.na(data) | (data == hdr$nodata_value)
        for (na.val in nodata.values)  na = na | (data==na.val)
        if (any(na)) data[na] = NA
    } else {
        data = matrix(NA,ncol=hdr$ncols,nrow=hdr$nrows)
        for (i in 1:hdr$nrows) {
            if (print == 2) cat(i, " ", ifelse(round(i/20)==i/20,"\n","") )
            if (print == 1) if (round(i/100)==i/100) cat(i, " ", ifelse(round(i/1000)==i/1000,"\n",""))
            x = scan(con, nlines = 1, quiet = TRUE, na.strings = na.strings)
            na = is.na(x) | (x == hdr$nodata_value)
            for (na.val in nodata.values)  na = na | (x==na.val)
            if (any(na)) x[na] = NA
            data[i,] = x
        }
    }
    if (print == 2) cat("\nDone!\n")
    if (print == 1) cat("\n")
    if (return.header) data = list( header = hdr, data = data )
    invisible(data)
}




#' @rdname read.ascii.grid
#' @name read.ascii.grid.header
#' @export
read.ascii.grid.header = function(file,...)
{
    if (is.character(file)) {
        file = default.file.extension(file,".asc")
        file = file(file,open="r")
        on.exit(close(file), add = TRUE)
    }
    hdr = scan(file, what=list(attribute="",value=numeric(0)), 
            nlines=6, quiet=TRUE, ...)
    hdr$attribute = tolower(hdr$attribute)
    res = hdr$value
    names(res) = hdr$attribute
    res = as.list(res)
    if (!is.null(res$xllcorner) & !is.null(res$yllcorner)) {
        res$xllcenter = res$xllcorner + res$cellsize / 2
        res$yllcenter = res$yllcorner + res$cellsize / 2
    } else if (!is.null(res$xllcenter) & !is.null(res$yllcenter)) {
        res$xllcorner = res$xllcenter - res$cellsize / 2
        res$yllcorner = res$yllcenter - res$cellsize / 2
    }
    return(res)
}

#' @rdname read.ascii.grid
#' @name read.sgrd
#' @export
read.sgrd = function( fname, return.header = TRUE, print = 0, 
                      nodata.values = c(), at.once = TRUE, prec = 7, ... )
{
    temp.fname = paste(tempfile(),".asc",sep="")
    res = rsaga.sgrd.to.esri( fname, temp.fname, prec=prec, format="ascii",
                              show.output.on.console=FALSE, intern=FALSE, ... )
    on.exit(unlink(temp.fname), add = TRUE)
    if (res==0) {
        data = read.ascii.grid( temp.fname, return.header=return.header,
                                print=print, nodata.values=nodata.values, at.once=at.once )
    } else
        stop("error converting the SAGA sgrd file to a temporary ASCII grid file")
    invisible(data)
}

#' @rdname read.ascii.grid
#' @name read.Rd.grid
#' @export
read.Rd.grid = function( fname, return.header = TRUE )
{
    fname = default.file.extension(fname,".Rd")
    load(fname)
    stopifnot(exists("data", envir=parent.frame()))
    if (is.list(data)) 
        stopifnot( (names(data)==c("header","data")) | (names(data)==c("data","header")) )
    if (return.header & !is.list(data)) { 
        warning("header missing")
        data = list(header=NA,data=data)
    } else if (!return.header & is.list(data)) 
        data = data$data
    invisible(data)
}

#' @rdname read.ascii.grid
#' @name write.ascii.grid
#' @export
write.ascii.grid = function( data, file, header = NULL, write.header = TRUE, 
                             digits, hdr.digits = 10, dec = ".", georef = "corner" ) 
{
    if (is.character(file)) {
        file = default.file.extension(file, ".asc")
        con = file(file,open="w")
        on.exit(close(con), add = TRUE)
    } else {
        if (!isOpen(file,"write"))
            stop("'file' must be a file name or a connection opened for writing")
    }
    if (is.list(data)) {
        stopifnot( ("data" %in% names(data)) )
        if (write.header & is.null(header)) {
            stopifnot("header" %in% names(data))
            header = data$header
        }
        data = data$data
    } else stopifnot(is.matrix(data))
    if (!missing(digits)) 
        data = round(data,digits=digits)
    if (write.header)  
        write.ascii.grid.header(con, header, dec=dec, georef=georef, hdr.digits=hdr.digits)
    utils::write.table(data, file=con, append=TRUE, quote=FALSE,
                na=as.character(header$nodata_value),
                row.names=FALSE, col.names=FALSE, dec=dec)
}

#' @rdname read.ascii.grid
#' @name write.ascii.grid.header
#' @export
write.ascii.grid.header = function(file, header, georef, dec=".", hdr.digits=10)
    # added digits argument - 2013-02-07
{
    if (missing(georef)) {
        # determine from the 'header' if georeferencing should refer
        # to corner or center of lower left grid cell:
        i.corner = min(c(Inf,grep("corner",tolower(names(header)))))
        i.center = min(c(Inf,grep("center",tolower(names(header)))))
        stopifnot(i.center!=i.corner) # this can only happen if header is corrupt
        georef = "corner"
        if (i.center < i.corner) georef = "center"
    } else {
        georef = match.arg(tolower(georef),choices=c("corner","center"))
    }
    # number of decimal places in header now determined by digits argument; 2013-02-07:
    my.fmt = paste("%-14s%-.",as.character(hdr.digits),"f",sep="")
    fmt = c("%-14s%-.0f", "%-14s%-.0f", my.fmt, my.fmt, 
            my.fmt, my.fmt)
    nm = c( "ncols", "nrows", paste(c("xll","yll"),georef,sep=""), "cellsize", "nodata_value" )
    if (is.character(file))  {
        file = default.file.extension(file,".asc")
        file = file(file, open="w")
        on.exit(close(file), add = TRUE)
    } else {
        if (!isOpen(file,"write"))
            stop("'file' must be a file name or a connection opened for reading")
    }
    for (i in 1:length(nm)) {
        entry = gsub(".", dec, sprintf(fmt[i],nm[i],as.numeric(header[[ nm[i] ]])), fixed=TRUE)
        write( entry, file=file, append=(i>1) )
    }
    invisible()
}

#' @rdname read.ascii.grid
#' @name write.sgrd
#' @export
write.sgrd = function( data, file, header = NULL, prec = 7,    
                       hdr.prec = 10, georef = "corner", ... )
    # 'georef' argument was missing - bug fixed 2008-05-02
    # hdr.prec argument added - 2013-02-07
{
    temp.fname = paste(tempfile(),".asc",sep="")
    write.ascii.grid( data = data, file = temp.fname, header = header, 
                      digits = prec, hdr.digits = hdr.prec, georef = georef )
    on.exit(unlink(temp.fname), add = TRUE)
    res = rsaga.esri.to.sgrd( in.grids = temp.fname, out.sgrds = file,
                              show.output.on.console = FALSE, intern = FALSE, ... )
    invisible(res)
}

#' @rdname read.ascii.grid
#' @name write.Rd.grid
#' @export
write.Rd.grid = function(data, file, header=NULL, write.header=TRUE, 
    compress=TRUE)
{
    file = default.file.extension(file,".Rd")
    if (is.list(data)) {
        stopifnot( ("data" %in% names(data)) )
        if (write.header & is.null(header)) {
            stopifnot("header" %in% names(data))
            header = data$header
        }
        data = data$data
    } else stopifnot(is.matrix(data))
    if (write.header)  data = list( header = header, data = data )
    save(data, file=file, ascii=FALSE, compress=compress)
}


#' Pick Variable from Spatial Dataset
#' 
#' These functions pick (i.e. interpolate without worrying too much about theory) values of a spatial variables from a data stored in a data.frame, a point shapefile, or an ASCII or SAGA grid, using nearest neighbor or kriging interpolation. \code{pick.from.points} and \code{[internal.]pick.from.ascii.grid} are the core functions that are called by the different wrappers.
#' @name pick.from.points
#' @param data data.frame giving the coordinates (in columns specified by \code{X.name, Y.name}) of point locations at which to interpolate the specified variables or grid values
#' @param src data.frame
#' @param shapefile point shapefile
#' @param pick variables to be picked (interpolated) from \code{src}; if missing, use all available variables, except those specified by \code{X.name} and \code{Y.name}
#' @param method interpolation method to be used; uses a partial match to the alternatives \code{"nearest.neighbor"} (currently the default) and \code{"krige"}
#' @param set.na logical: if a column with a name specified in \code{pick} already exists in \code{data}, how should it be dealt with? \code{set.na=FALSE} (default) only overwrites existing data if the interpolator yields a non-\code{NA} result; \code{set.na=TRUE} passes \code{NA} values returned by the interpolator on to the results data.frame
#' @param radius numeric value specifying the radius of the local neighborhood to be used for interpolation; defaults to 200 map units (presumably meters), or, in the functions for grid files, \code{2.5*cellsize}.
#' @param nmin numeric, for \code{method="krige"} only: see \code{\link[gstat]{krige}} function in package \pkg{gstat}
#' @param nmax numeric, for \code{method="krige"} only: see \code{\link[gstat]{krige}} function in package \pkg{gstat}
#' @param sill numeric, for \code{method="krige"} only: the overall sill parameter to be used for the variogram
#' @param range numeric, for \code{method="krige"} only: the variogram range
#' @param nugget numeric, for \code{method="krige"} only: the nugget effect
#' @param model for \code{method="krige"} only: the variogram model to be used for interpolation; defaults to a spherical variogram with parameters specified by the \code{range}, \code{sill}, and \code{nugget} arguments; see \code{\link[gstat]{vgm}} in package \pkg{gstat} for details
#' @param log logical vector, specifying for each variable in \code{pick} if interpolation should take place on the logarithmic scale (default: \code{FALSE})
#' @param X.name name of the variable containing the x coordinates
#' @param Y.name name of the variable containing the y coordinates
#' @param cbind logical: shoud the new variables be added to the input data.frame (\code{cbind=TRUE}, the default), or should they be returned as a separate vector or data.frame? \code{cbind=FALSE}
#' @param file file name (relative to \code{path}, default file extension \code{.asc}) of an ASCII grid from which to pick a variable, or an open connection to such a file
#' @param path optional path to \code{file}
#' @param varname character string: a variable name for the variable interpolated from grid file \code{file} in \code{pick.from.*.grid}; if missing, variable name will be determined from \code{file}name by a call to \code{\link{create.variable.name}}
#' @param prefix an optional prefix to be added to the \code{varname}
#' @param nodata.values numeric vector specifying grid values that should be converted to \code{NA}; in addition to the values specified here, the nodata value given in the input grid's header will be used
#' @param at.once logical: should the grid be read as a whole or line by line? \code{at.once=FALSE} is useful for processing large grids that do not fit into memory; the argument is currently by default \code{FALSE} for \code{method="nearest.neighbour"}, and it currently MUST be \code{TRUE} for all other methods (in these cases, \code{TRUE} is the default value); piecewise processing with \code{at.once=FALSE} is always faster than processing the whole grid \code{at.once}
#' @param quiet logical: provide information on the progress of grid processing on screen? (only relevant if \code{at.once=FALSE} and \code{method="nearest.neighbour"})
#' @param nlines numeric: stop after processing \code{nlines} lines of the input grid; useful for testing purposes
#' @param filename character: name of a SAGA grid file, default extension \code{.sgrd}
#' @param prec numeric, specifying the number of digits to be used in converting a SAGA grid to an ASCII grid in \code{pick.from.saga.grid}
#' @param na.strings passed on to \code{\link{scan}}
#' @param env list: RSAGA geoprocessing environment created by \code{\link{rsaga.env}}
#' @param show.output.on.console a logical (default: \code{FALSE}), indicates whether to capture the output of the command and show it on the R console (see \code{\link{system}}, \code{\link{rsaga.geoprocessor}}).
#' @param nsplit split the data.frame \code{data} in \code{nsplit} disjoint subsets in order to increase efficiency by using \code{\link[plyr]{ddply}} in package \pkg{plyr}. The default seems to perform well in many situations.
#' @param parallel logical (default: \code{FALSE}): enable parallel processing; requires additional packages such as \pkg{doSNOW} or \pkg{doMC}. See example below and \code{\link[plyr]{ddply}}
#' @param ... arguments to be passed to \code{pick.from.points}, and to \code{internal.pick.from.ascii.grid} in the case of \code{pick.from.ascii.grid}
#' @details \code{pick.from.points} interpolates the variables defined by \code{pick} in the \code{src} data.frame to the locations provided by the \code{data} data.frame. Only nearest neighbour and ordinary kriging interpolation are currently available. This function is intended for 'data-rich' situations in which not much thought needs to be put into a geostatistical analysis of the spatial structure of a variable. In particular, this function is supposed to provide a simple, 'quick-and-dirty' interface for situations where the \code{src} data points are very densely distributed compared to the \code{data} locations.
#'
#' \code{pick.from.shapefile} is a front-end of \code{pick.from.points} for point shapefiles.
#' 
#' \code{pick.from.ascii.grid} retrieves data values from an ASCII raster file using either nearest neighbour or ordinary kriging interpolation. The latter may not be possible for large raster data sets because the entire grid needs to be read into an R matrix. Split-apply-combine strategies are used to improve efficiency and allow for parallelization.
#' 
#' The optional parallelization of \code{pick.from.ascii.grid} computation requires the use of a \emph{parallel backend} package such as \pkg{doSNOW} or \pkg{doMC}, and the parallel backend needs to be registered before calling this function with \code{parallel=TRUE}. The example section provides an example using \pkg{doSNOW} on Windows. I have seen 25-40% reduction in processing time by parallelization in some examples that I ran on a dual core Windows computer.
#' 
#' \code{pick.from.ascii.grids} performs multiple \code{pick.from.ascii.grid} calls. File \code{path} and \code{prefix} arguments may be specific to each \code{file} (i.e. each may be a character vector), but all interpolation settings will be the same for each \code{file}, limiting the flexibility a bit compared to individual \code{pick.from.ascii.grid} calls by the user. \code{pick.from.ascii.grids} currently processes the files sequentially (i.e. parallelization is limited to the \code{pick.from.ascii.grid} calls within this function).
#' 
#' \code{pick.from.saga.grid} is the equivalent to \code{pick.from.ascii.grid} for SAGA grid files. It simply converts the SAGA grid \code{file} to a (temporary) ASCII raster file and applies \code{pick.from.ascii.grid}.
#' 
#' \code{internal.pick.from.ascii.grid} is an internal 'workhorse' function that by itself would be very inefficient for large data sets \code{data}. This function is called by \code{pick.from.ascii.grid}, which uses a split-apply-combine strategy implemented in the \pkg{plyr} package.
#' 
#' @return If \code{cbind=TRUE}, columns with the new, interpolated variables are added to the input data.frame \code{data}.
#'
#' If \code{cbind=FALSE}, a data.frame only containing the new variables is returned (possibly coerced to a vector if only one variable is processed).
#' 
#' @references Brenning, A. (2008): Statistical geocomputing combining R and SAGA:  The example of landslide susceptibility analysis with generalized additive models. In: J. Boehner, T. Blaschke, L. Montanarella (eds.), SAGA - Seconds Out (= Hamburger Beitraege zur Physischen Geographie und Landschaftsoekologie, 19), 23-32.
#'
#' @author Alexander Brenning
#' @note \code{method="krige"} requires the \pkg{gstat} package.
#' 
#' \code{pick.from.shapefile} requires the \pkg{shapefiles} package.
#' 
#' The nearest neighbour interpolation currently randomly breaks ties if \code{pick.from.points} is used, and in a deterministic fashion (rounding towards greater grid indices, i.e. toward south and east) in the grid functions.
#'
#' @seealso  \code{\link{grid.to.xyz}}, %\code{\link{vgm}}, \code{\link{krige}}, \code{\link{read.ascii.grid}}, \code{\link{write.ascii.grid}}
#' @examples 
#' \dontrun{
#' # assume that 'dem' is an ASCII grid and d a data.frame with variables x and y
#' pick.from.ascii.grid(d, "dem")
#' # parallel processing on Windows using the doSNOW package:
#' require(doSNOW)
#' registerDoSNOW(cl <- makeCluster(2, type = "SOCK")) # DualCore processor
#' pick.from.ascii.grid(d, "dem", parallel = TRUE)
#' # produces two (ignorable) warning messages when using doSNOW
#' # typically 25-40% faster than the above on my DualCore notebook
#' stopCluster(cl)
#' }
#' 
#' \dontrun{
#' # use the meuse data for some tests:
#' require(gstat)
#' data(meuse)
#' data(meuse.grid)
#' meuse.nn = pick.from.points(data=meuse.grid, src=meuse, 
#'     pick=c("cadmium","copper","elev"), method="nearest.neighbour")
#' meuse.kr = pick.from.points(data=meuse.grid, src=meuse,
#'     pick=c("cadmium","copper","elev"), method="krige", radius=100)
#' # it does make a difference:
#' plot(meuse.kr$cadmium,meuse.nn$cadmium)
#' plot(meuse.kr$copper,meuse.nn$copper)
#' plot(meuse.kr$elev,meuse.nn$elev)
#' }
#' @keywords spatial
#' @export
pick.from.points = function(data, src, pick, 
    method = c("nearest.neighbour","krige"), set.na = FALSE,
    radius = 200, nmin = 0, nmax = 100,
    sill = 1, range = radius, nugget = 0,
    model = vgm(sill-nugget,"Sph",range=range,nugget=nugget),
    log = rep(FALSE,length(pick)),
    X.name = "x", Y.name = "y", cbind = TRUE )
{
    method = match.arg(method)
    if (!is.data.frame(data))  data = as.data.frame(data)
    if (missing(pick)) {
        pick = colnames(src)
        pick = pick[ (pick != X.name) & (pick != Y.name) ]
    }
    nc = rep(NA,length(pick))
    for (p in 1:length(pick)) {
        colnm = colnames(data)
        if (any(colnm==pick[p])) {
            nc[p] = which(colnm==pick[p])
            if (set.na)  data[,nc[p]] = rep(NA,nrow(data))
        } else {
            data = cbind(data,rep(NA,nrow(data)))
            colnames(data) = c(colnm,pick[p])
            nc[p] = ncol(data)
        }
    }
    stopifnot(!any(is.na(nc)))
    
    # prepare the source data.frame:
    src = src[ !is.na(src[,X.name]) & !is.na(src[,Y.name]) , ]
    if (nrow(src)==0) {
        if (!cbind) data = data[,pick]
        return(data)
    }
    the.src = src
    rm(src)
    
    if (method=="krige") {
        loc = stats::as.formula(paste("~",X.name,"+",Y.name))
        for (p in 1:length(pick)) {
            form = stats::as.formula(paste(pick[p],"~ 1"))
            src = the.src[ !is.na(the.src[,pick[p]]) , ]
            if (nrow(src)==0) next
            krg = krige(
                form, locations=loc, data=src, newdata=data,
                model=model, nmax=nmax, nmin=nmin,
                maxdist=radius ) $ var1.pred
            sel = !is.na(krg)
            data[sel,pick[p]] = krg[sel]
        }
    } else if (method=="nearest.neighbour") {
        for (i in 1:nrow(data)) {
            if (is.na(data[i,X.name]) | is.na(data[i,Y.name])) next
            di = sqrt((the.src[,X.name]-data[i,X.name])^2+(the.src[,Y.name]-data[i,Y.name])^2)
            mindi = min(di)
            if ((mindi) > radius) next
            
            wh = which(di == mindi)
            if (length(wh)>1) wh = sample(wh)[1]
            for (p in 1:length(pick))
                data[i,pick[p]] = the.src[wh,pick[p]]
        }
    }
    if (!cbind) data = data[,pick]
    return(data)
}


#' @rdname pick.from.points
#' @name pick.from.shapefile
#' @export
pick.from.shapefile = function(data, shapefile, X.name="x", Y.name="y", ...)
{
    shapefile = set.file.extension(shapefile,"")
    shapefile = substr(shapefile,1,nchar(shapefile)-1) # remove "." at the end
    src = read.shapefile(shapefile)
    src = add.xy(src)
    src = src$dbf[[1]]
    if (X.name != "XCOORD") {
        if (X.name %in% names(src)) {
            src[,X.name] = src[,"XCOORD"]
            src = src[,names(src)!="XCOORD"]
        }
    }
    if (Y.name != "YCOORD") {
        if (Y.name %in% names(src)) {
            src[,Y.name] = src[,"YCOORD"]
            src = src[,names(src)!="YCOORD"]
        }
    }
    data = pick.from.points(data,src,X.name=X.name,Y.name=Y.name,...)
    return(data)
}

#' @rdname pick.from.points
#' @name pick.from.ascii.grid
#' @export
pick.from.ascii.grid = function( data, file, path = NULL, varname = NULL, prefix = NULL,
                                 method = c("nearest.neighbour","krige"), cbind = TRUE,
                                 parallel = FALSE, nsplit, quiet = TRUE, ... )
{
    method = match.arg(method)
    
    # TO DO: parallel implementation not currently working:
    # it won't find the 'file' file unless the full path is specified
    parallel = FALSE
    
    if (missing(nsplit)) {
        if (method == "krige") {
            nsplit = 1 + parallel
        } else {
            nsplit = ceiling(nrow(data) / 1500)
            if (parallel)  nsplit = max(2, nsplit)
        }
    }
    
    if (nsplit == 1) {
        return( internal.pick.from.ascii.grid(data = data, file = file, path = path, varname = varname, 
                                              prefix = prefix, method = method, quiet = quiet, ...))
    } else {
        progress = "none"
        if (parallel)  quiet = TRUE
        if (nrow(data) >= 1000 & !quiet) {
            progress = ifelse(Sys.info()["sysname"]=="Windows" & .Platform$GUI == "Rgui", "win", "text")
            quiet = TRUE
        }
        PICKSPLIT = floor( seq(0, nsplit-0.001, length = nrow(data)) )
        if (cbind) {
            op = options(warn=-1)
            on.exit(options(op))
            data = ddply( data, .variables = .(PICKSPLIT), .fun = internal.pick.from.ascii.grid,
                          file = file, path = path, varname = varname, prefix = prefix, method = method, 
                          quiet = quiet, cbind = cbind, ...,
                          .progress = progress, .parallel = parallel )
            options(op)
            data$PICKSPLIT = NULL
            return(data)
        } else {
            op = options(warn=-1)
            on.exit(options(op))
            res = dlply( data, .variables = .(PICKSPLIT), .fun = internal.pick.from.ascii.grid,
                         file = file, path = path, varname = varname, prefix = prefix, method = method, 
                         quiet = quiet, cbind = cbind, ...,
                         .progress = progress, .parallel = parallel )
            options(op)
            ###print(str(res))
            res = unlist(res, use.names = FALSE)
            ###print(str(res))
            return(res)
        }
    }
}

#' @rdname pick.from.points
#' @name pick.from.ascii.grids
#' @export
pick.from.ascii.grids = function( data, file, path = NULL, varname = NULL, prefix = NULL,
                                  cbind = TRUE, quiet = TRUE, ... )
{
    if (!is.null(path)) {
        if (length(path) == 1) path = rep(path, length(file))
        stopifnot(length(path) == length(file))
    }
    if (!is.null(varname)) {
        stopifnot(length(varname) == length(file))
    }
    if (!is.null(prefix)) {
        if (length(prefix) == 1) prefix = rep(prefix, length(file))
        stopifnot(length(prefix) == length(file))
    }
    
    if (length(file) == 1) {
        return( pick.from.ascii.grid( data = data, file = file, path = path, varname = varname, 
                                      prefix = prefix, cbind = cbind, quiet = quiet, ...) )
    }
    
    if (is.null(varname)) {
        if (is.character(file)) {
            varname = unname( sapply(file, RSAGA::create.variable.name) )
        } else {
            if (cbind) {
                stop("'varname' must be specified unless 'file' is a character string with the filename")
            } else varname = paste("X", c(1:length(file)), sep = "")
        }
    }
    
    # add a prefix to the variable names?
    if (!is.null(prefix))
        for (i in 1:length(file))
            if (prefix!="")
                varname[i] = paste(prefix[i],varname[i],sep=".")
            
            for (i in 1:length(file)) {
                if (!quiet) cat("Processing file '", file[i], "' (", i, " of ", length(file), ")...\n", sep="")
                res = pick.from.ascii.grid( data = data, file = file[i], path = path[i], varname = varname[i],
                                            prefix = prefix[i], cbind = cbind, quiet = TRUE, ...)
                if (cbind) {
                    data = res
                } else {
                    if (i == 1) {
                        RES = res
                    } else RES = cbind(RES, res)
                }
            }
            if (!cbind) {
                data = as.data.frame(RES)
                colnames(data) = varname
            }
            return(data)
}

#' @rdname pick.from.points
#' @name internal.pick.from.ascii.grid
#' @export
internal.pick.from.ascii.grid = function( data, file, 
    path = NULL, varname = NULL, prefix = NULL,
    method = c("nearest.neighbour","krige"),
    nodata.values = c(-9999,-99999), at.once, quiet = TRUE,
    X.name = "x", Y.name = "y", nlines = Inf,
    cbind = TRUE, range, radius, na.strings = "NA", ... )
{
    stopifnot(is.data.frame(data))
    stopifnot( X.name %in% colnames(data) )
    stopifnot( Y.name %in% colnames(data) )

    # determine variable name from file name if 'varname' is missing:
    if (is.null(varname)) {
        if (is.character(file)) {
            varname = RSAGA::create.variable.name(file)
        } else {
            if (cbind) {
                stop("'varname' must be specified unless 'file' is a character string with the filename")
            } else
                varname = paste("TEMP",X.name,Y.name,sep=".")
        }
    }

    # add a prefix to the variable?
    if (!is.null(prefix)) if (prefix!="")
        varname = paste(prefix,varname,sep=".")

    method = match.arg(method)
    
    if (missing(at.once))
        at.once = (method != "nearest.neighbour")
    
    if (is.character(file)) {
        file = RSAGA::default.file.extension(file,".asc")
        if (!is.null(path)) if (path!="") file = file.path(path,file)
        if (!file.exists(file)) stop("file ",file," not found")
        con = file(file,open="r")
        on.exit(close(con), add = TRUE)
    } else {
        con = file
        if (!isOpen(con,"read"))
            stop("'file' must be a file name, or a connection opened for reading")
    }

    # prepare 'data': add new column if necessary
    the.colnames = colnames(data)
    if (varname %in% the.colnames) {
        nc = which(the.colnames==varname)
    } else {
        data = cbind( data, rep(NA,nrow(data) ) )
        colnames(data) = c( the.colnames, varname )
        nc = ncol(data)
    }

    if (method=="krige")
    {
        if (!at.once)
            warning("row-by-row processing of grids is not yet implemented for kriging interpolation\n",
                "trying to process the whole grid at once...")
        src = RSAGA::read.ascii.grid(con, nodata.values = nodata.values, 
                na.strings = na.strings)
        src = RSAGA::grid.to.xyz(src, colnames=c(X.name,Y.name,varname))
        if (missing(radius)) radius = 2.5 * hdr$cellsize
        if (missing(range)) range = radius
        if (range > radius) radius = range
        data = RSAGA::pick.from.points(data, src, pick = varname,
            X.name=X.name, Y.name = Y.name,
            method="krige", range = range, radius = radius,...)

    } else if (method=="nearest.neighbour")
    {
        hdr = RSAGA::read.ascii.grid.header(con)
        nodata.values = unique(c(nodata.values,hdr$nodata_value))

        select = cbind( 1 + round( (data[,X.name] - (hdr$xllcorner+hdr$cellsize/2)) / hdr$cellsize ),
                        1 + round( hdr$nrows - (data[,Y.name] - (hdr$yllcorner-hdr$cellsize/2)) / hdr$cellsize )  )
    
        if (any(!is.na(select)))
        {
            nr = nrow(select)
            nlines = max( 1, min( hdr$nrows, max(select[,2],na.rm=TRUE), nlines ) )
        
            if (!at.once)
            {
                for (i in 1:nlines) {
                    if (!quiet & ((i/10)==floor(i/10))) {
                        cat(i," ")
                        if ( (i/100) == floor(i/100) )  cat("\n")
                    }
                    v = scan(con, nlines = 1, quiet = TRUE, na.strings = na.strings)
                    if (length(v) < hdr$ncols) {
                        warning("grid row too short - corrupt grid file? filling with NA's...")
                        v = c(v, rep(nodata.values[1], hdr$ncols-length(v)))
                    } else if (length(v) > hdr$ncols) {
                        warning("grid row too long - corrupt grid file? ignoring redundant data...")
                        v = v[ 1:hdr$ncols ]
                    }
                    ass = ( select[,2] == i )
                    if (any(ass)) {
                        for (na in nodata.values) v[ v==na ] = NA
                        for (j in which(ass)) {
                            if ((select[j,1]>=1) & (select[j,1]<=hdr$ncols))
                                data[j,nc] = v[ select[j,1] ]
                        }
                        ###if (!quiet) cat("\n matches:",which(ass),"\n")
                    }
                }
            } else # if (at.once)
            {
                v = utils::read.table(con, na.strings = na.strings)
                for (na in nodata.values) v[ v==na ] = NA
                for (j in 1:nr) {
                    if ( all(select[j,]>=1) & all(select[j,]<=c(ncol(v),nrow(v))) )
                        data[j,nc] = v[select[j,2],select[j,1]]
                }
            }
        } else {
            warning("all 'data' points are outside grid area")
        }
    } # end if (method=="nearest.neighbour")
    
    if (!cbind) data = data[,nc]
##print(str(data)); cat("----\n")
    return(data)
}


#' @rdname pick.from.points
#' @name pick.from.saga.grid
#' @export
pick.from.saga.grid = function( data, filename, path, varname, 
                                prec = 7, show.output.on.console = FALSE, env = rsaga.env(), ... )
{
    if (!missing(path)) if (path!="") filename = file.path(path,filename)
    temp.asc = paste(tempfile(),".asc",sep="")
    if (missing(varname)) varname = create.variable.name(filename)
    rsaga.sgrd.to.esri(filename, temp.asc, format = "ascii",
                       prec = prec, show.output.on.console = show.output.on.console,
                       env = env)
    on.exit(unlink(temp.asc), add = TRUE)
    data = pick.from.ascii.grid(data, temp.asc, varname = varname, ...)
    invisible(data)
}


#' Convert Grid Matrix to (x,y,z) data.frame
#' 
#' Convert a grid matrix to a (x,y,z) data.frame.
#' @name grid.to.xyz
#' @param data grid data: either a grid data matrix, or a list with components \code{data} (a matrix with the grid data) and \code{header} (the grid header information); see \code{\link{read.ascii.grid}} for details 
#' @param header optional list giving grid header information; see \code{\link{read.ascii.grid}} for details
#' @param varname character: name to be assigned to the column with the z values in the output data.frame
#' @param colnames names to be given to the columns corresponding to the x and y coordinates and the grid variable in the output data.frame
#' @return a data.frame with three columns (names are specified in the \code{colnames} argument) giving the x and y coordinates and the attribute values at the locations given by the grid \code{data}.
#' @seealso \code{\link{read.ascii.grid}}, \code{\link{pick.from.ascii.grid}}
#' @examples 
#' \dontrun{
#' d = read.ascii.grid("dem")
#' xyz = grid.to.xyz(d,varname="elevation")
#' str(xyz)
#' }
#' @keywords spatial
#' @export
grid.to.xyz = function(data,header,varname="z",colnames=c("x","y",varname)) {
    if (missing(header)) {
        if (is.list(data)) {
            header = data$header
        } else {
            header = list(
                ncols = ncol(data),
                nrows = nrow(data),
                xllcenter = 0,
                yllcenter = 0,
                cellsize = 1,
                xllcorner = -0.5,
                yllcorner = -0.5 )
        }
    }
    if (is.list(data)) data = data$data
    data = data.frame(
        x = header$xllcenter + rep( seq( 0, (header$ncols-1) * header$cellsize, by=header$cellsize ), header$nrows ),
        y = rep( header$yllcenter + seq( (header$nrows-1) * header$cellsize, 0, by=-header$cellsize ), each=header$ncols ),
        z = as.vector(t(data)) )
    colnames(data) = colnames
    invisible(data)
}




#' Pick Center Value from Matrix
#'
#' Pick the value in the center of a square matrix. Auxiliary function to be used by functions called by \code{\link{focal.function}}.
#' @name centervalue
#' @param x a square matrix
#' @details See for example the code of \code{\link{resid.median}}.
#' @seealso \code{\link{focal.function}}, \code{\link{resid.median}}
#' @examples 
#' ( m <- matrix( round(runif(9,1,10)), ncol=3 ) )
#' centervalue(m)
#' @keywords utilities
#' @export
centervalue = function(x) {
    i = ceiling(ncol(x) / 2)
    return(x[i,i])
}



#' Residual Median and Quantile Filters for Grids
#'
#' These functions use the median and other quantiles to describe the difference between a grid value and its neighborhood. They are designed for use with \code{\link{focal.function}}.
#' @name resid.median
#' @param x a square matrix with the grid data from the moving window, possibly containing \code{NA} values
#' @param probs numeric vector of probabilities in [0,1] to be passed to \code{\link{quantile}}
#' @details These functions are designed for being called by \code{\link{focal.function}}, which repeatedly passes the contents of a square or circular moving window to these functions.
#' 
#' The \code{resid.median} function rests the value of the central grid cell from the median of the whole moving window. Thus, in terms of topography, a positive residual median indicates that this grid cell stands out compared to its surroundings. \code{resid.quantile} gives more flexibility in designing such residual attributes.
#' @return If \code{x} is provided, a numeric vector of length 1 (\code{resid.median}), 3 (\code{resid.minmedmax} and \code{resid.quartiles}), or \code{length(probs)} (\code{resid.quantile}).
#'
#' If \code{x} is missing, a character vector of same length giving suggested variable (or file) names, such as \code{"rmed"}. See \code{\link{focal.function}} for details.
#' @seealso \code{\link{focal.function}}, \code{\link{quantile}}, \code{\link{median}}, \code{\link{centervalue}}
#' @keywords spatial
#' @export
resid.median = function(x) {
    if (missing(x)) return("rmed")
    return( stats::median(x,na.rm=TRUE) - centervalue(x) )
}

#' @rdname resid.median
#' @name resid.minmedmax
#' @export
resid.minmedmax = function(x) {
    if (missing(x)) return(c("rmin","rmed","rmax"))
    return( c(min(x,na.rm=TRUE),stats::median(x,na.rm=TRUE),max(x,na.rm=TRUE)) - centervalue(x) )
}

#' Relative Topographic Position
#' 
#' \code{relative.position} and \code{relative.rank} are used with \code{\link{focal.function}} to determine the relative value of a grid cell compared to its surroundings, either on a metric scale or based on ranks.
#' @name relative.position
#' @param x a square matrix with the grid data from the moving window, possibly containing \code{NA} values
#' @param ties.method see \code{\link{rank}}
#' @return If \code{x} is provided, a numeric value in the interval [0,1] is returned.
#'
#' If \code{x} is missing, a character vector of same length giving suggested variable (or file) names, here \code{"relpos"} and \code{"relrank"}, respectively. See \code{\link{focal.function}} for details.
#' @seealso \code{\link{focal.function}}, \code{\link{rank}}, \code{\link{centervalue}} 
#' @examples 
#' m = matrix( round(runif(9,1,10)), ncol=3 )
#' print(m)
#' relative.position(m)
#' relative.rank(m)
#' \dontrun{
#' focal.function("dem",fun=relative.rank,radius=5)
#' focal.function("dem",fun=relative.position,radius=5)
#' relrank = as.vector(read.ascii.grid("relrank")$data)
#' relpos  = as.vector(read.ascii.grid("relpos")$data)
#' plot(relpos,relrank,pch=".")
#' cor(relpos,relrank,use="complete.obs",method="pearson")
#' }
#' @keywords spatial
#' @export
relative.position = function(x) {
    if (missing(x)) return("relpos")
    zmin = min(x,na.rm=TRUE)
    zmax = max(x,na.rm=TRUE)
    return( (centervalue(x) - zmin) / (zmax - zmin) )
}

#' @rdname resid.median
#' @name resid.quantile
#' @export
resid.quantile = function(x,probs) {
    if (missing(x)) return(NULL)
    return(stats::quantile(x-centervalue(x),probs=probs,na.rm=TRUE,names=FALSE))
}

#' @rdname resid.median
#' @name resid.quartiles
#' @export
resid.quartiles = function(x) {
    if (missing(x)) return(c("r25","r50","r75"))
    return(stats::quantile(x-centervalue(x),probs=c(0.25,0.5,0.75),na.rm=TRUE,names=FALSE))
}

#' @rdname relative.position
#' @name relative.rank
#' @export
relative.rank = function(x,ties.method="average") {
    if (missing(x)) return("relrank")
    x = as.vector(x)
    n = sum(!is.na(x))
    return( (rank(x,ties.method=ties.method)[(length(x)+1)/2]-1) / (n-1))
}


#' Wind Shelter Index
#' 
#' \code{wind.shelter} is a function to be used with \code{\link{focal.function}} to calculate a topographic wind shelter index from a digital elevation model, which is a proxy for snow accumulation on the lee side of topographic obstacles. \code{wind.shelter.prep} performs some preparatory calculations to speed up repeated calls to \code{wind.shelter}.
#' @name wind.shelter
#' @param x square matrix of elevation data
#' @param prob numeric: quantile of slope values to be used in computing the wind shelter index; if \code{NULL}, use \code{max} (equivalent to \code{prob=1})
#' @param control required argument: the result of a call to \code{wind.shelter.prep}
#' @param radius radius (>1) of circle segment to be used (number of grid cells, not necessarily an integer)
#' @param direction wind direction: direction from which the wind originates; North = 0 = \code{2*pi}, clockwise angles.
#' @param tolerance directional tolerance
#' @param cellsize grid cellsize
#' @details \code{wind.shelter} implements a wind shelter index used by Plattner et al. (2004) for modeling snow accumulation patterns on a glacier in the Austrian Alps. It is a modified version of the algorithm of Winstral et al. (2002). The wind shelter index of Plattner et al. (2004) is defined as:
#'
#'    \code{Shelter index(S) = arctan( max( (z(x0)-z(x)) / |x0-x| : x in S ) ),}
#'    
#' where \code{S = S(x0,a,da,d)} is the set of grid nodes within a distance \code{<=d} from \code{x0}, only considering grid nodes in directions between \code{a-da} and \code{a+da} from \code{x0}.
#'
#' The present implementation generalizes this index by replacing \code{max} by the \code{quantile} function; the \code{max} function is used if \code{prob=NULL}, and the same result is obtained for \code{prob=1} using the \code{quantile} function.
#' @return The function \code{wind.shelter} returns the wind shelter index as described above if a numeric matrix \code{x} is provided. If it is missing, it returns the character string \code{"windshelter"}.
#'  
#' \code{wind.shelter.prep} returns a list with components \code{mask} and \code{dist}. Both are square matrices with \code{2*(ceiling(radius)+1)} columns and rows:
#'   \item{mask}{indicates which grid cell in the moving window is within the specified circle segment (value \code{FALSE}) or not (\code{TRUE})}
#'   \item{dist}{the precomputed distances of a grid cell to the center of the moving window, in map units}
#' @references Plattner, C., Braun, L.N., Brenning, A. (2004): Spatial variability of snow accumulation on Vernagtferner, Austrian Alps, in winter 2003/2004. Zeitschrift fuer Gletscherkunde und Glazialgeologie, 39: 43-57.
#'
#' Winstral, A., Elder, K., Davis, R.E. (2002): Spatial snow modeling of wind-redistributed snow using terrain-based parameters. Journal of Hydrometeorology, 3: 524-538.
#' @author Alexander Brenning
#' @note The wind shelter index only makes sense if elevation is measured in the same units as the horizontal map units used for the \code{cellsize} argument (i.e. usually meters).
#'
#' \code{wind.shelter} and \code{wind.shelter.prep} do not restrict the calculation to a circular area; this is done by \code{\link{focal.function}} when used in combination with that function (assuming \code{search.mode="circle"}).
#'
#' Note that the present definition of the wind shelter index returns negative values for surfaces that are completely exposed toward the specified direction. This may make sense if interpreted as a "wind exposure index", or it might be appropriate to set negative wind shelter values to 0.
#' @seealso  \code{\link{focal.function}}, \code{\link{quantile}} 
#' @examples 
#' # Settings used by Plattner et al. (2004):
#' ctrl = wind.shelter.prep(6,-pi/4,pi/12,10)
#' \dontrun{focal.function("dem.asc",fun=wind.shelter,control=ctrl,
#'     radius=6,search.mode="circle")}
#' @keywords spatial
#' @export
wind.shelter = function(x,prob=NULL,control) {
    if (missing(x)) return("windshelter")
    if (missing(control)) stop("need 'control' argument - call 'wind.shelter.prep' first")
    ctr = centervalue(x)
    x[control$mask] = NA
    res = NA
    if (!all(is.na(x))) {
        x = atan((x-ctr)/control$dist)
        if (is.null(prob)) {
            res = max(x,na.rm=TRUE)
        } else res = stats::quantile(x,probs=prob,na.rm=TRUE)
    }
    return(res)
}




#' @rdname wind.shelter
#' @name wind.shelter.prep
#' @export
wind.shelter.prep = function(radius,direction,tolerance,cellsize=90) {
    nc = nr = 2*ceiling(radius)+1
    mask = matrix(TRUE,ncol=nc,nrow=nr)
    for (j in 1:nc) {
        for (i in 1:nr) {
            if ((i==j) & (i==((nr+1)/2))) next
            xy = c( j-(nc+1)/2, (nr+1)/2-i )
            xy = xy / sqrt(xy[1]^2+xy[2]^2)
            if ( xy[2]>0)  a = asin(xy[1])  else a = pi - asin(xy[1])
            if (a < 0) a = a + 2*pi
            d = abs(direction-a)
            if (d>2*pi) d = d-2*pi
            d = min(d,2*pi-d)
            if (d<=tolerance) mask[i,j] = FALSE
        }
    }
    dist = matrix(NA,ncol=nc,nrow=nr)
    for (i in 1:nr) for (j in 1:nc) {
        xy = c( j-(nc+1)/2, (nr+1)/2-i )
        dist[i,j] = sqrt(xy[1]^2+xy[2]^2) * cellsize
    }
    list( mask = mask, dist = dist )
}




#' Local and Focal Grid Functions
#' 
#' \code{focal.function} cuts out square or circular moving windows from a grid (matrix) and applies a user-defined matrix function to calculate e.g. a terrain attribute or filter the grid. The function is suitable for large grid files as it can process them row by row. \code{local.function} represents the special case of a moving window of radius 1. Users can define their own functions operating on moving windows, or use simple functions such as \code{median} to define filters.
#' @name focal.function
#' @param in.grid file name of input ASCII grid, relative to \code{in.path}
#' @param in.factor.grid optional file name giving a gridded categorical variables defining zones; zone boundaries are used as breaklines for the moving window (see Details)
#' @param out.grid.prefix character string (optional), defining a file name prefix to be used for the output file names; a dash (\code{-}) will separate the prefix and the \code{varnames} 
#' @param path path in which to look for \code{in.grid} and write output grid files; see also \code{in.path} and \code{out.path}, which overwrite \code{path} if they are specified
#' @param in.path path in which to look for \code{in.grid} (defaults to \code{path})
#' @param out.path path in which to write output grid files; defaults to \code{path}
#' @param fun a function, or name of a function, to be applied on the moving window; see Details
#' @param varnames character vector specifying the names of the variable(s) returned by \code{fun}; if missing, \code{focal.function} will try to determine the varnames from \code{fun} itself, or from a call to \code{fun} if this is a function (see Details)
#' @param radius numeric value specifying the (circular or square) radius  of the moving window; see \code{is.pixel.radius} and \code{search.mode}; note that all data within distance \code{<=radius} will be included in the moving window, not \code{<radius}.
#' @param is.pixel.radius logical: if \code{TRUE} (default), the \code{radius} will be interpreted as a (possibly non-integer) number of pixels; if \code{FALSE}, it is interpreted as a radius measured in the grid (map) units.
#' @param valid.range numeric vector of length 2, specifying minimum and maximum valid values read from input file; all values \code{<valid.range[1]} or \code{>valid.range[1]} will be converted to \code{NA}.
#' @param nodata.values numeric vector: any values from the input grid file that should be converted to \code{NA}, in addition to the nodata value specified in the grid header
#' @param out.nodata.value numeric: value used for storing \code{NA}s in the output file(s); if missing, use the same nodata value as specified in the header of the input grid file
#' @param na.strings passed on to \code{\link{scan}}
#' @param search.mode character, either \code{"circle"} (default) for a circular search window, or \code{"square"} for a squared one.
#' @param digits numeric, specifying the number of digits to be used for output grid file.
#' @param hdr.digits numeric, specifying the number of digits to be used for the header of the output grid file (default: 10; see \code{\link{write.ascii.grid.header}}).
#' @param dec character, specifying the decimal mark to be used for input and output.
#' @param quiet If \code{TRUE}, gives some output (\code{"*"}) after every 10th line of the grid file and when the job is done.
#' @param nlines Number of lines to be processed; useful for testing purposes.
#' @param mw.to.vector logical: Should the content of the moving window be coerced (from a matrix) to a vector?
#' @param mw.na.rm logical: Should \code{NA}s be removed from moving window prior to passing the data to \code{fun}? Only applicable when \code{mw.to.vector=TRUE}.
#' @param \dots Arguments to be passed to \code{fun}; \code{local.function}: arguments to be passed to  \code{focal.function}.
#' @details \code{focal.function} passes a square matrix of size \code{2*radius+1} to the function \code{fun} if \code{mw.to.vector=FALSE} (default), or a vector of length \code{<=(2*radius+1)^2} if \code{mw.to.vector=TRUE}. This matrix or vector will contain the content of the moving window, which may possibly contain \code{NA}s even if the \code{in.grid} has no nodata values, e.g. due to edge effects. If \code{search.mode="circle"}, values more than \code{radius} units (pixels or grid units, depending on \code{is.pixel.radius}) away from the center pixel / matrix entry will be set to \code{NA}. In addition, \code{valid.range}, \code{nodata.values}, and the nodata values specified in the \code{in.grid} are checked to assign further \code{NA}s to pixels in the moving window. Finally, if \code{in.factor.grid} specifies zones, all pixels in the moving window that belong to a different zone than the center pixel are set to \code{NA}, or, in other words, zone boundaries are used as breaklines.
#'
#' The function \code{fun} should return a single numeric value or a numeric vector. As an example, the function \code{\link{resid.minmedmax}} returns the minimum, median and maximum of the difference between the values in the moving window and the value in the center grid cell. In addition to the (first) argument receiving the moving window data, \code{fun} may have additional arguments; the \code{...} argument of \code{focal.function} is passed on to \code{fun}. \code{\link{resid.quantile}} is a function that uses this feature.
#'
#' Optionally, \code{fun} should support the following feature: If no argument is passed to it, then it should return a character vector giving variable names to be used for naming the output grids. The call \code{\link{resid.minmedmax}()}, for example, returns \code{c("rmin","rmed","rmax")}; this vector must have the same length as the numeric vector returned when moving window data is passed to the function. This feature is only used if no \code{varnames} argument is provided. Note that the result is currently being \code{\link{abbreviate}}d to a length of 6 characters.
#'
#' Input and output file names are built according to the following schemes:
#'
#' Input:   \code{[<in.path>/]<in.grid>}
#'
#' Zones:   \code{[<in.path>/]<in.factor.grid>} (if specified)
#'
#' Output:  \code{[<out.path>/][<out.grid.prefix>-]<varnames>.asc}
#'
#' For the input files, \code{.asc} is used as the default file extension, if it is not specified by the user.
#' 
#' @return \code{focal.function} and \code{local.function} return the character vector of output file names.
#' @references Brenning, A. (2008): Statistical geocomputing combining R and SAGA: The example of landslide susceptibility analysis with generalized additive models.  In: J. Boehner, T. Blaschke, L. Montanarella (eds.), SAGA - Seconds Out (= Hamburger Beitraege zur Physischen Geographie und Landschaftsoekologie, 19), 23-32.
#' @author Alexander Brenning
#' @note These functions are not very efficient ways of calculating e.g. (focal) terrain attributes compared to for example the SAGA modules, but the idea is that you can easily specify your own functions without starting to mess around with C code. For example try implementing a median filter as a SAGA module... or just use the code shown in the example!
#' @seealso \code{\link{multi.focal.function}}, \code{\link{multi.local.function}}, \code{\link{resid.median}}, \code{\link{resid.minmedmax}}, \code{\link{relative.position}}, \code{\link{resid.quantile}}, \code{\link{resid.quartiles}}, \code{\link{relative.rank}},  \code{\link{wind.shelter}}, \code{\link{create.variable.name}}
#' @examples 
#' \dontrun{
#' # A simple median filter applied to dem.asc:
#' gapply("dem","median",radius=3)
#' # Same:
#' #focal.function("dem",fun="median",radius=3,mw.to.vector=TRUE,mw.na.rm=TRUE)
#' # See how the filter has changed the elevation data:
#' d1 = as.vector(read.ascii.grid("dem")$data)
#' d2 = as.vector(read.ascii.grid("median")$data)
#' hist(d1-d2,br=50)
#' }
#' # Wind shelter index used by Plattner et al. (2004):
#' \dontrun{
#' ctrl = wind.shelter.prep(6,-pi/4,pi/12,10)
#' focal.function("dem",fun=wind.shelter,control=ctrl,
#'     radius=6,search.mode="circle")
#' }
#' # Or how about this, if "aspect" is local terrain exposure:
#' \dontrun{
#' gapply("aspect","cos") # how "northerly-exposed" is a pixel?
#' gapply("aspect","sin") # how "easterly-exposed" is a pixel?
#' # Same result, but faster:
#' focal.function("aspect",fun=function(x) c(cos(x),sin(x)), varnames=c("cos","sin"))
#' }
#' @keywords spatial
#' @export
focal.function = function( in.grid, in.factor.grid, out.grid.prefix,
    path=NULL, in.path=path, out.path=path,
    fun, varnames,
    radius=0, is.pixel.radius=TRUE,
    na.strings = "NA",
    valid.range=c(-Inf,Inf), nodata.values=c(), out.nodata.value, 
    search.mode=c("circle","square"),
    digits=4, hdr.digits=10, dec=".", quiet=TRUE, nlines=Inf,
    mw.to.vector = FALSE, mw.na.rm = FALSE, ... )
{
    if (radius > 0) {
        search.mode = match.arg(search.mode)
        if (mw.na.rm & !mw.to.vector)
            warning("'mw.na.rm=TRUE' only meaningful if moving window matrix is\n",
                "converted to a vector ('mw.to.vector=TRUE')")
    }
    
    # prepare input file:
    if (!is.null(in.path)) if (in.path!="")
        in.grid = file.path(in.path,in.grid)
    in.grid = default.file.extension(in.grid,".asc")
    in.file = file(in.grid,open="r")
    on.exit(close(in.file), add = TRUE)
    in.hdr = read.ascii.grid.header(in.file,dec=dec)
    nodata.values = unique(c(nodata.values,in.hdr$nodata_value))
    nlines = max( 1, min( c(nlines,in.hdr$nrows), na.rm=TRUE ) )

    if (missing(in.factor.grid)) in.factor.grid = NULL
    if ((radius<=0) & !is.null(in.factor.grid)) {
        warning("'in.factor.grid' is ignored - only meaningful for 'radius>0'")
        in.factor.grid = NULL
    }
    if (!is.null(in.factor.grid)) {
        in.factor.grid = file.path(in.path,in.factor.grid)
        in.factor.grid = default.file.extension(in.factor.grid,".asc")
        in.factor.file = file(in.factor.grid,open="r")
        on.exit(close(in.factor.file), add = TRUE)
        in.factor.hdr = read.ascii.grid.header(in.factor.file,dec=dec)
        if (in.hdr$ncols != in.factor.hdr$ncols |
            in.hdr$nrows != in.factor.hdr$nrows |
            in.hdr$cellsize != in.factor.hdr$cellsize |
            in.hdr$xllcorner != in.factor.hdr$xllcorner |
            in.hdr$yllcorner != in.factor.hdr$yllcorner)
            stop("input grid and factor grid must have same extent and cellsize")
    }

    # build output filenames:
    if (missing(varnames)) {
        # check if the function will return a vector with variable names
        # when called without arguments:
        varnames = try(do.call(fun,list()),silent=TRUE)
        if (missing(varnames) || class(varnames) == "try-error") {
            if (is.character(fun)) {
                varnames = gsub(".","",fun,fixed=TRUE)
            } else if (is.function(fun)) {
                varnames = deparse(substitute(fun))
            } else stop("unable to determine 'varnames' from 'fun'")
            varnames = abbreviate(varnames,6)
        }
    }
    if (missing(out.grid.prefix)) out.grid.prefix = ""
    if (is.null(out.grid.prefix)) out.grid.prefix = ""
    stopifnot(length(varnames) == length(unique(varnames)))
    do.paste = (varnames!="") & (out.grid.prefix!="")
    out.filenames = paste( out.grid.prefix, c("","_")[do.paste+1], varnames, sep="" )
    out.filenames = default.file.extension(out.filenames,".asc")
    if (!is.null(out.path)) if (out.path!="")
        out.filenames = file.path(out.path,out.filenames)
    if (any(out.filenames==in.grid)) stop("one of the output file names is identical to the input file name")

    # prepare output files:
    N.out = length(out.filenames)
    out.files = as.list(1:N.out)
    out.hdr = in.hdr
    if (missing(out.nodata.value)) out.nodata.value = in.hdr$nodata_value
    out.hdr$nodata_value = out.nodata.value
    for (k in 1:N.out) {
        out.files[[k]] = file(out.filenames[k],open="w")
        write.ascii.grid.header(out.files[[k]],out.hdr,dec=dec,hdr.digits=hdr.digits)
    }
    on.exit( for (k in 1:N.out) close(out.files[[k]]), add=TRUE )
    fmt = paste("%.",digits,"f",sep="")

    if (radius <= 0) {
        # Apply 'fun' as a local function:
    
        # Process one line at a time:
        for (i in 1:nlines) {
            if (!quiet) if ((i %% 10)==0) cat("*")
            if (!quiet) if ((i %% 100)==0) cat("\n")
            
            # Read one line at a time:
            v0 = scan(in.file, nlines = 1, quiet = TRUE, dec = dec, 
                    na.strings = na.strings)
            if (length(v0) != in.hdr$ncols) {
                warning("grid line does not have NCOLS values")
                v0 = c( v0, rep(NA,in.hdr$ncols-length(v0)) )
            }
            for (na in nodata.values)  v0[ v0==na ] = NA
            v0[ v0 < valid.range[1] ] = NA
            v0[ v0 > valid.range[2] ] = NA
             
#            # With plyr package instead of for loop:   
#            require(plyr)
#            mycall = function(x,...) do.call(fun,list(x,...))
#            res = t(laply(.data = as.list(v0), .fun = mycall, .drop = FALSE, .parallel = parallel, ...))
#            # ...but for some reason it doesn't work with .parallel = TRUE

            res = matrix(NA,ncol=in.hdr$ncols,nrow=N.out)
            for (j in 1:in.hdr$ncol) {
                r = do.call(fun,list(v0[j],...))
                res[,j] = r
            }

            res[ is.na(res) ] = out.nodata.value
            for (k in 1:N.out) {
                txt = paste(sprintf(fmt,res[k,]),collapse=" ")
                if (dec!=".") txt = gsub(".",dec,txt,fixed=TRUE)
                writeLines(txt,con=out.files[[k]])
            }
        }
    
    } else { # if (radius > 0)
    
        if (!is.pixel.radius) radius = radius / in.hdr$cellsize
        exact.radius = radius
        radius = ceiling(radius)
    
        # 'v' is a matrix that will receive a set of rows copied from the grid;
        # it must be a bit wider than the grid so the moving window can move over
        # it without having to worry about edge effects:
        v = matrix( NA, ncol=in.hdr$ncols+2*radius, nrow=2*radius+1 )
        # 'fac': same for in.factor.grid, if available:
        if (!is.null(in.factor.grid))
            fac = matrix( NA, ncol=in.hdr$ncols+2*radius, nrow=2*radius+1 )
        # 'f' will be the mask of a moving window in case of a circular window:    
        if (search.mode=="circle") {
            f = matrix(FALSE,ncol=2*radius+1,nrow=2*radius+1)
            for (i in ((-1)*radius):radius)
                for (j in ((-1)*radius):radius)
                    if (sqrt(i^2+j^2) > exact.radius)
                        f[ i+radius+1, j+radius+1 ] = TRUE
        }
        
        # the look-ahead step:
        for (i in (radius+1):(2*radius)) {
            v[i+1,] = c( rep(NA,radius), 
                scan(in.file, nlines = 1, quiet = TRUE, dec = dec,
                    na.strings = na.strings), 
                rep(NA,radius) )
            if (!is.null(in.factor.grid))
                fac[i+1,] = c( rep(NA,radius), 
                    scan(in.factor.file, nlines = 1, quiet = TRUE,
                        dec = dec, na.strings = na.strings), 
                    rep(NA,radius) )
        }
        # Process nodata values:
        for (na in nodata.values)  v[ v==na ] = NA
        v[ v < valid.range[1] ] = NA
        v[ v > valid.range[2] ] = NA
        # Process nodata values of the factor grid:
        if (!is.null(in.factor.grid)) {
            fac[ fac==in.factor.hdr$nodata_value ] = NA
            v[ is.na(fac) ] = NA
        }
        
        # Process the grid line by line:
        for (i in 1:nlines) {
            if (!quiet) if ((i %% 10)==0) cat("*")
            if (!quiet) if ((i %% 100)==0) cat("\n")
            
            if (i <= nlines - radius) {
                # Read a line from the grid file:
                v0 = scan(in.file, nlines = 1, quiet = TRUE, dec = dec,
                        na.strings = na.strings)
                if (length(v0) != in.hdr$ncols) { # check if corrupt
                    warning("grid line does not have NCOLS values")
                    v0 = c( v0, rep(NA,ncol(v)-length(v0)) )
                }
                # process all the nodata values:
                for (na in nodata.values)  v0[ v0==na ] = NA
                v0[ v0 < valid.range[1] ] = NA
                v0[ v0 > valid.range[2] ] = NA

                # Read a line from the factor grid:
                if (!is.null(in.factor.grid)) {
                    fac0 = scan(in.factor.file, nlines = 1, quiet = TRUE,
                        na.strings = na.strings)
                    if (length(fac0) != in.factor.hdr$ncols) {
                        warning("factor grid line does not have NCOLS values")
                        fac0 = c( fac0, rep(NA,ncol(fac)-length(fac0)) )
                    }
                    fac0[ fac0 == in.factor.hdr$nodata_value ] = NA
                    # Pass NA's on to the grid itself:
                    v0[ is.na(fac0) ] = NA
                }
            } else {
                v0 = rep(NA,in.hdr$ncols)
                if (!is.null(in.factor.grid))  fac0 = v0
            }
            
            # Add new line to the look-ahead buffer:
            v = rbind( v[2:(2*radius+1),], t(c( rep(NA,radius), v0, rep(NA,radius) )) )
            if (!is.null(in.factor.grid))
                fac = rbind( fac[2:(2*radius+1),], t(c( rep(NA,radius), fac0, rep(NA,radius) )) )
            
            # Apply the 'fun'ction to each grid column:
            res = matrix(NA,ncol=in.hdr$ncol,nrow=N.out)
            for (j in 1:in.hdr$ncol) {
                w = v[,j:(j+2*radius)]
                if (search.mode=="circle")  w[f] = NA
                if (!is.null(in.factor.grid)) {
                    facw = fac[,j:(j+2*radius)]
                    the.fac = centervalue(facw)
                    if (is.na(the.fac)) {
                        w = NA
                    } else
                        w[ facw != the.fac ] = NA
                }
                if (!all(is.na(w))) {
                    if (mw.to.vector) {
                        w = as.vector(w)
                        if (mw.na.rm) w = w[!is.na(w)]
                    }
                    r = do.call(fun,list(w,...))
                    res[,j] = r
                }
            }
            res[ is.na(res) ] = out.nodata.value
            for (k in 1:N.out) {
                txt = paste(sprintf(fmt,res[k,]),collapse=" ")
                if (dec!=".") txt = gsub(".",dec,txt,fixed=TRUE)
                writeLines(txt,con=out.files[[k]])
            }
        }
    } # end if (radius > 0)
    
    if (!quiet)  cat("\nDone.\n")
    return(out.filenames)
}




#' @rdname focal.function
#' @name gapply
#' @export
gapply = function(in.grid,fun,varnames,mw.to.vector=TRUE,mw.na.rm=TRUE,...) {
    # build output filenames:
    if (missing(varnames)) {
        # check if the function will return a vector with variable names
        # when called without arguments:
        varnames = try(do.call(fun,list()),silent=TRUE)
        if (class(varnames) == "try-error") {
            if (is.character(fun)) {
                varnames = gsub(".","",fun,fixed=TRUE)
            } else if (is.function(fun)) {
                varnames = deparse(substitute(fun))
            } else stop("unable to determine 'varnames' from 'fun'")
            varnames = abbreviate(varnames,6)
        }
    }
    focal.function(in.grid=in.grid,fun=fun,varnames=varnames,
        mw.to.vector=mw.to.vector,mw.na.rm=mw.na.rm,...)
}

#' @rdname focal.function
#' @name local.function
#' @export
local.function = function( ... ) {
    focal.function(..., radius=0,
        in.factor.grid=NULL, search.mode=NULL, is.pixel.radius=NULL,
        mw.to.vector=FALSE, mw.na.rm=FALSE )
}





#' Local and Focal Grid Function with Multiple Grids as Inputs
#'
#' \code{multi.focal.function} cuts out square or circular moving windows from a stack of grids (matrices) and applies a user-defined matrix function that takes multiple arguments to this data. \code{multi.local.function} is a more efficiently coded special case of moving windows of size 0, i.e. functions applied to individual grid cells of a stack of grids. This is especially useful for applying \code{predict} methods of statistical models to a stack of grids containing the explanatory variables (see Examples and \code{\link{grid.predict}}). The function is suitable for large grid files as it can process them row by row; but it may be slow because one call to the focal function is generated for each grid cell.
#' 
#' @name multi.focal.function
#' @param in.grids character vector: file names of input ASCII grids, relative to \code{in.path}; \code{in.grid.prefix} will be used as a prefix to the file name if specified; default file extension: \code{.asc}
#' @param in.factor.grid optional file name giving a gridded categorical variables defining zones; zone boundaries are used as breaklines for the moving window (see Details)
#' @param in.grid.prefix character string (optional), defining a file name prefix to be used for the input file names; a dash (\code{-}) will separate the prefix and the \code{in.varnames}
#' @param out.grid.prefix character string (optional), defining a file name prefix to be used for the output file names; a dash (\code{-}) will separate the prefix and the \code{out.varnames}
#' @param path path in which to look for \code{in.grids} and write output grid files; see also \code{in.path} and \code{out.path}, which overwrite \code{path} if they are specified
#' @param in.path path in which to look for \code{in.grids} (defaults to \code{path})
#' @param out.path path in which to write output grid files; defaults to \code{path}
#' @param fun a function, or name of a function, to be applied on the moving window; see Details; \code{fun} is expected to accept named arguments with the names given by \code{in.varnames}; \code{\link{grid.predict}} is a wrapper function that can be used for applying a model's \code{predict} method to a stack of grids; see Details. In \code{multi.local.function}, \code{fun} must be able to process  arguments that are vectors of equal length (e.g., a vector of 50 slope angles, another vector of 50 elevation values, etc.).
#' @param in.varnames character vector: names of the variables corresponding to the \code{in.grids}; if missing, same as \code{in.grids}; if specified, must have the same length and order as \code{in.grids}
#' @param out.varnames character vector specifying the name(s) of the variable(s) returned by \code{fun}; if missing, \code{multi.focal.function} will try to determine the varnames from \code{fun} itself, or or from a call to \code{fun} if this is a function (see Details)
#' @param radius numeric value specifying the (circular or square) radius  of the moving window; see \code{is.pixel.radius} and \code{search.mode}; note that all data within distance \code{<=radius} will be included in the moving window, not \code{<radius}.
#' @param is.pixel.radius logical: if \code{TRUE} (default), the \code{radius} will be interpreted as a (possibly non-integer) number of pixels; if \code{FALSE}, it is interpreted as a radius measured in the grid (map) units.
#' @param valid.ranges optional list of length \code{length(in.grids)} with numeric vector of length 2, specifying minimum and maximum valid values read from input file; all values \code{<valid.ranges[[i]][1]} or \code{>valid.ranges[[i]][1]} will be converted to \code{NA}.
#' @param nodata.values numeric vector: any values from the input grid file that should be converted to \code{NA}, in addition to the nodata value specified in the grid header
#' @param out.nodata.value numeric: value used for storing \code{NA}s in the output file(s); if missing, use the same nodata value as specified in the header of the input grid file
#' @param search.mode character, either \code{"circle"} (default) for a circular search window, or \code{"square"} for a squared one.
#' @param digits numeric, specifying the number of digits to be used for output grid file.
#' @param hdr.digits numeric, specifying the number of digits to be used for the header of the output grid file (default: 10; see \code{\link{write.ascii.grid.header}}).
#' @param dec character, specifying the decimal mark to be used for input and output.
#' @param quiet If \code{FALSE}, gives some output (\code{"*"}) after every 10th line of the grid file and when the job is done.
#' @param nlines Number of lines to be processed; useful for testing purposes.
#' @param na.action function: determines if/how \code{NA} values are omitted from the stack of input variables; use \code{\link{na.exclude}} (default) or \code{\link{na.pass}} if \code{fun} can handle \code{NA} values correctly
#' @param mw.to.vector logical: Should the content of the moving window be coerced (from a matrix) to a vector?
#' @param mw.na.rm logical: Should \code{NA}s be removed from moving window prior to passing the data to \code{fun}? Only applicable when \code{mw.to.vector=TRUE}.
#' @param pass.location logical: Should the x,y coordinates of grid points (center of grid cells) be passed to \code{fun}? If \code{TRUE}, two additional arguments named arguments \code{x} and \code{y} are passed to \code{fun}; NOTE: This currently only works for \code{radius=0}, otherwise a warning is produced and \code{pass.location} is reset to \code{FALSE}.
#' @param na.strings passed on to \code{\link{scan}}
#' @param \dots Arguments to be passed to \code{fun}; \code{local.function}: arguments to be passed to  \code{focal.function}.
#' @details \code{multi.local.function} is probably most useful for applying the \code{predict} method of a fitted model to a grids representing the predictor variables. An example is given below and in more detail in Brenning (2008) (who used \code{multi.focal.function} for the same purpose); see also \code{\link{grid.predict}}.
#'
#' \code{multi.local.function} is essentially the same as \code{multi.focal.function} for \code{radius=0}, but coded MUCH more efficiently. (The relevant code will eventually migrate into \code{multi.focal.function} as well, but requires further testing.) Applying a GAM to the data set of Brenning (2008) takes about 1/100th the time with \code{multi.local.function} compared to \code{multi.focal.function}.
#' 
#' \code{multi.focal.function} extends \code{\link{focal.function}} by allowing multiple input grids to be passed to the focal function \code{fun} operating on moving windows. It passes square matrices of size \code{2*radius+1} to the function \code{fun} if \code{mw.to.vector=FALSE} (default), or a vector of length \code{<=(2*radius+1)^2} if \code{mw.to.vector=TRUE}; one such matrix or vector per input grid will be passed to \code{fun} as an argument whose name is specified by \code{in.varnames}.
#' 
#' These matrices or vectors will contain the content of the moving window, which may possibly contain \code{NA}s even if the \code{in.grid} has no nodata values, e.g. due to edge effects. If \code{search.mode="circle"}, values more than \code{radius} units (pixels or grid units, depending on \code{is.pixel.radius}) away from the center pixel / matrix entry will be set to \code{NA}. In addition, \code{valid.range}, \code{nodata.values}, and the nodata values specified in the \code{in.grid} are checked to assign further \code{NA}s to pixels in the moving window. Finally, if \code{in.factor.grid} specifies zones, all pixels in the moving window that belong to a different zone than the center pixel are set to \code{NA}, or, in other words, zone boundaries are used as breaklines.
#'
#' The function \code{fun} should return a single numeric value or a numeric vector, such as a regression result or a vector of class probabilities returned by a soft classifier. In addition to the named arguments receiving the moving window data, \code{fun} may have additional arguments; the \code{...} argument of \code{focal.function} is passed on to \code{fun}. \code{\link{grid.predict}} uses this feature.
#'
#' Optionally, \code{fun} should support the following feature: If no argument is passed to it, then it should return a character vector giving variable names to be used for naming the output grids.
#'
#' For the input files, \code{.asc} is used as the default file extension, if it is not specified by the user.
#'
#' See \code{\link{focal.function}} for details.
#'
#' @return \code{multi.focal.function} returns the character vector of output file names.
#' @references Brenning, A. (2008): Statistical geocomputing combining R and SAGA: The example of landslide susceptibility analysis with generalized additive models. In: J. Boehner, T. Blaschke, L. Montanarella (eds.), SAGA - Seconds Out (= Hamburger Beitraege zur Physischen Geographie und Landschaftsoekologie, 19), 23-32.
#' @author Alexander Brenning
#' @note \code{multi.focal.function} can do all the things \code{\link{focal.function}} can do.
#' @seealso \code{\link{focal.function}}, \code{\link{grid.predict}}
#' @examples 
#' \dontrun{
#' # Assume that d is a data.frame with point observations
#' # of a numerical response variable y and predictor variables
#' # a, b, and c.
#' # Fit a generalized additive model to y,a,b,c.
#' # We want to model b and c as nonlinear terms:
#' require(gam)
#' fit <- gam(y ~ a + s(b) + s(c), data = d)
#' multi.local.function(in.grids = c("a", "b", "c"),
#'     out.varnames = "pred",
#'     fun = grid.predict, fit = fit )
#'     # Note that the 'grid.predict' uses by default the
#'     # predict method of 'fit'.
#' # Model predictions are written to a file named pred.asc
#' }
#'
#' \dontrun{
#' # A fake example of a logistic additive model:
#' require(gam)
#' fit <- gam(cl ~ a + s(b) + s(c), data = d, family = binomial)
#' multi.local.function(in.grids = c("a", "b", "c"),
#'     out.varnames = "pred",
#'     fun = grid.predict, fit = fit,
#'     control.predict = list(type = "response") )
#'     # 'control.predict' is passed on to 'grid.predict', which
#'     # dumps its contents into the arguments for 'fit''s
#'     # 'predict' method.
#' # Model predictions are written to a file named pred.asc
#' }
#' @keywords spatial
#' @export
multi.focal.function = function( 
    in.grids, in.grid.prefix, in.factor.grid, 
    out.grid.prefix,
    path = NULL, in.path = path, out.path = path,
    fun, in.varnames, out.varnames,
    radius = 0, is.pixel.radius = TRUE,
    na.strings = "NA",
    valid.ranges, nodata.values = c(), out.nodata.value, 
    search.mode = c("circle","square"),
    digits = 4, hdr.digits = 10, dec = ".", quiet = TRUE, nlines = Inf,
    mw.to.vector = FALSE, mw.na.rm = FALSE, pass.location = FALSE, 
    ... )
{
    if (radius > 0) {
        search.mode = match.arg(search.mode)
        if (mw.na.rm & !mw.to.vector)
            warning("'mw.na.rm=TRUE' only meaningful if moving window matrix is\n",
                "converted to a vector ('mw.to.vector=TRUE')")
    }
    
    # build input filenames:
    if (missing(in.grid.prefix)) in.grid.prefix = ""
    if (is.null(in.grid.prefix)) in.grid.prefix = ""
    if (missing(in.varnames)) {
        in.varnames = in.grids
    } else if (missing(in.grids)) {
        in.grids = in.varnames
    }
    stopifnot(length(in.varnames) == length(unique(in.varnames)))
    stopifnot(length(in.grids) == length(unique(in.grids)))
    stopifnot(length(in.varnames) == length(in.grids))
    do.paste.in = (in.varnames!="") & (in.grid.prefix!="")
    in.filenames = paste( in.grid.prefix, c("","_")[do.paste.in+1], in.grids, sep="" )
    in.filenames = default.file.extension(in.filenames,".asc")
    if (!is.null(in.path)) if (any(in.path != ""))
        in.filenames = file.path(in.path, in.filenames)

    # prepare input files:
    N.in = length(in.filenames)
    in.files = in.hdrs = nodata.vals = as.list(1:N.in)
    for (k in 1:N.in) {
        in.files[[k]] = file(in.filenames[k],open="r")
        in.hdrs[[k]] = read.ascii.grid.header(in.files[[k]],dec=dec)
        nodata.vals[[k]] = unique(c(nodata.values,in.hdrs[[k]]$nodata_value))
        if (k > 1) {
            if ( in.hdrs[[k]]$cellsize != in.hdrs[[1]]$cellsize |
                  in.hdrs[[k]]$ncols != in.hdrs[[1]]$ncols |
                  in.hdrs[[k]]$nrows != in.hdrs[[1]]$nrows )
                stop("incompatible input grids")
        }
    }
    on.exit( for (k in 1:N.in) close(in.files[[k]]) ) # add = TRUE
    in.hdr = in.hdrs[[1]]
    nlines = max( 1, min( c(nlines,in.hdr$nrows), na.rm=TRUE ) )

    if (missing(in.factor.grid)) in.factor.grid = NULL
    if ((radius<=0) & !is.null(in.factor.grid)) {
        warning("'in.factor.grid' is ignored - only meaningful for 'radius>0'")
        in.factor.grid = NULL
    }
    if (!is.null(in.factor.grid)) {
        in.factor.grid = file.path(in.path,in.factor.grid)
        in.factor.grid = default.file.extension(in.factor.grid,".asc")
        in.factor.file = file(in.factor.grid,open="r")
        on.exit(close(in.factor.file),add=TRUE)
        in.factor.hdr = read.ascii.grid.header(in.factor.file,dec=dec)
        if (in.hdr$ncols != in.factor.hdr$ncols |
            in.hdr$nrows != in.factor.hdr$nrows |
            in.hdr$cellsize != in.factor.hdr$cellsize )
            stop("input grid and factor grid must have same extent and cellsize")
    }

    # build output filenames:
    if (missing(out.varnames)) {
        # check if the function will return a vector with variable names
        # when called without arguments:
        out.varnames = try(do.call(fun,list()),silent=TRUE)
        if (missing(out.varnames) || class(out.varnames) == "try-error") {
            if (is.character(fun)) {
                out.varnames = gsub(".","",fun,fixed=TRUE)
            } else if (is.function(fun)) {
                out.varnames = deparse(substitute(fun))
            } else stop("unable to determine 'out.varnames' from 'fun'")
            out.varnames = abbreviate(out.varnames,6)
        }
    }
    if (missing(out.grid.prefix)) out.grid.prefix = ""
    if (is.null(out.grid.prefix)) out.grid.prefix = ""
    stopifnot(length(out.varnames) == length(unique(out.varnames)))
    do.paste = (out.varnames!="") & (out.grid.prefix!="")
    out.filenames = paste( out.grid.prefix, c("","_")[do.paste+1], 
                    out.varnames, sep="" )
    out.filenames = default.file.extension(out.filenames,".asc")
    if (!is.null(out.path)) if (out.path!="")
        out.filenames = file.path(out.path,out.filenames)
    if (any(out.filenames %in% in.filenames))
        stop("one of the output file names is equal to an input file name")

    # prepare output files:
    N.out = length(out.filenames)
    out.files = as.list(1:N.out)
    out.hdr = in.hdr
    if (missing(out.nodata.value)) out.nodata.value = in.hdr$nodata_value
    out.hdr$nodata_value = out.nodata.value
    for (k in 1:N.out) {
        out.files[[k]] = file(out.filenames[k],open="w")
        write.ascii.grid.header(out.files[[k]],out.hdr,dec=dec,hdr.digits=hdr.digits)
    }
    on.exit( for (k in 1:N.out) close(out.files[[k]]), add=TRUE )
    
    if (missing(valid.ranges)) {
        valid.ranges = list()
        for (k in 1:N.in) valid.ranges[[k]] = c(-Inf, Inf)
    }

    fmt = paste("%.",digits,"f",sep="")
    loc = NULL

    if (radius <= 0) {
        # Apply 'fun' as a local function:
    
        # Process one line at a time:
        for (i in 1:nlines) {
            if (!quiet) if ((i %% 10)==0) cat("*")
            if (!quiet) if ((i %% 100)==0) cat("\n")
            
            y.coord = in.hdr$yllcenter + (in.hdr$nrows - i) * in.hdr$cellsize
            
            # Read one line at a time, file by file:
            vl0 = as.list(1:N.in)
            for (k in 1:N.in) {
                vl0[[k]] = scan(in.files[[k]], nlines = 1, quiet = TRUE, 
                    dec = dec, na.strings = na.strings)
                if (length(vl0[[k]]) != in.hdr$ncols) {
                    warning("grid line does not have NCOLS values")
                    vl0[[k]] = c( vl0[[k]], 
                            rep(NA, in.hdr$ncols - length(vl0[[k]])) )
                }
                for (na in nodata.vals[[k]]) 
                    vl0[[k]][ vl0[[k]] == na ] = NA
                vl0[[k]][ vl0[[k]] < valid.ranges[[k]][1] ] = NA
                vl0[[k]][ vl0[[k]] > valid.ranges[[k]][2] ] = NA
            }
                            
            res = matrix(NA, ncol = in.hdr$ncols, nrow = N.out)

            for (j in 1:in.hdr$ncol) {
                # Pass the (x,y) coordinates to the function?
                if (pass.location) {
                    x.coord = in.hdr$xllcenter + (j-1) * in.hdr$cellsize
                    loc = list( location = c(x = x.coord, y = y.coord) )
                }
                
                args = as.list(1:N.in)
                skip = FALSE
                for (k in 1:N.in)
                    skip = skip | is.na(args[[k]] <- vl0[[k]][j])
                if (!skip) {
                    names(args) = in.varnames
                    args = c( args, loc, alist(...) )
                    r = do.call(fun, args)
                    res[,j] = r
                }
            }
            res[ is.na(res) ] = out.nodata.value
            for (k in 1:N.out) {
                txt = paste( sprintf(fmt, res[k,]), collapse = " " )
                if (dec != ".") txt = gsub(".", dec, txt, fixed = TRUE)
                writeLines(txt,con = out.files[[k]])
            }
        }
    
    } else { # if (radius > 0)
    
        if (pass.location) {
            pass.location = FALSE
            warning("'pass.location=TRUE' is currently only implemented for 'radius=0'\n")
            # to do: set up moving window matrices with x and y coordinates, respectively??
        }
    
        if (!is.pixel.radius) radius = radius / in.hdr$cellsize
        exact.radius = radius
        radius = ceiling(radius)
    
        # 'vl' is a list of matrices, each of which
        # will receive a set of rows copied from the grid;
        # it must be a bit wider than the grid so the moving window can move over
        # it without having to worry about edge effects:
        vl = list(1:N.in)
        for (k in 1:N.in)
            vl[[k]] = matrix( NA, ncol = in.hdr$ncols + 2*radius, 
                                  nrow = 2*radius + 1 )
        # 'fac': same for in.factor.grid, if available:
        if (!is.null(in.factor.grid))
            fac = matrix( NA, ncol = in.hdr$ncols + 2*radius, 
                              nrow = 2*radius + 1 )
        # 'f' will be the mask of a moving window in case of a circular window:    
        if (search.mode=="circle") {
            f = matrix(FALSE, ncol = 2*radius + 1, nrow = 2*radius + 1)
            for (i in (-radius):radius)
                for (j in (-radius):radius)
                    if (sqrt(i^2+j^2) > exact.radius)
                        f[ i + radius + 1, j + radius + 1 ] = TRUE
        }
        
        # the look-ahead step:
        for (k in 1:N.in) {
            for (i in (radius+1):(2*radius)) {
                vl[[k]][i+1,] = c( rep(NA, radius), 
                        scan(in.files[[k]], nlines = 1, quiet = TRUE, dec = dec,
                            na.strings = na.strings), 
                        rep(NA, radius) )
                if (k == 1) {
                    if (!is.null(in.factor.grid)) {
                        fac[i+1,] = c( rep(NA, radius), 
                            scan(in.factor.file, nlines = 1, quiet = TRUE,
                                na.strings = na.strings), 
                            rep(NA, radius) )
                    }
                }
            }
            # Process nodata values:
            for (na in nodata.vals[[k]])  vl[[k]][ vl[[k]]==na ] = NA
            vl[[k]][ vl[[k]] < valid.ranges[[k]][1] ] = NA
            vl[[k]][ vl[[k]] > valid.ranges[[k]][2] ] = NA
            # Process nodata values of the factor grid:
            if (!is.null(in.factor.grid)) {
                if (k == 1)
                    fac[ fac == in.factor.hdr$nodata_value ] = NA
                vl[[k]][ is.na(fac) ] = NA
            }
        }
        
        # Process the grid line by line:
        for (i in 1:nlines) {
            if (!quiet) if ((i %% 10)==0) cat("*")
            if (!quiet) if ((i %% 100)==0) cat("\n")
            
            y.coord = in.hdr$yllcenter + (in.hdr$nrows - i) * in.hdr$cellsize

            vl0 = as.list(1:N.in)

            if (i <= nlines - radius) {
                # Read a line from the grid file:
                for (k in 1:N.in) {
                    vl0[[k]] = scan(in.files[[k]], nlines = 1, quiet = TRUE, 
                            dec = dec, na.strings = na.strings)
                    if (length(vl0[[k]]) != in.hdr$ncols) { # check if corrupt
                        warning("grid line does not have NCOLS values")
                        vl0[[k]] = c( vl0[[k]], rep(NA, in.hdr$ncols - length(vl0[[k]])) )
                    }
                    # process all the nodata values:
                    for (na in nodata.vals[[k]])  vl0[[k]][ vl0[[k]]==na ] = NA
                    vl0[[k]][ vl0[[k]] < valid.ranges[[k]][1] ] = NA
                    vl0[[k]][ vl0[[k]] > valid.ranges[[k]][2] ] = NA
                }

                # Read a line from the factor grid:
                if (!is.null(in.factor.grid)) {
                    fac0 = scan(in.factor.file, nlines = 1, quiet = TRUE,
                        na.strings = na.strings)
                    if (length(fac0) != in.factor.hdr$ncols) {
                        warning("factor grid line does not have NCOLS values")
                        fac0 = c( fac0, rep(NA, in.hdr$ncols - length(fac0)) )
                    }
                    fac0[ fac0 == in.factor.hdr$nodata_value ] = NA
                    # Pass NA's on to the grid itself:
                    for (k in 1:N.in)
                        vl0[[k]][ is.na(fac0) ] = NA
                }
            } else {
                for (k in 1:N.in)
                    vl0[[k]] = rep(NA, in.hdr$ncols)
                if (!is.null(in.factor.grid))  fac0 = vl0
            }

            # Add new line to the look-ahead buffer:
            for (k in 1:N.in)
                vl[[k]] = rbind( vl[[k]][2:(2*radius+1),], t(c( rep(NA,radius), vl0[[k]], rep(NA,radius) )) )
            if (!is.null(in.factor.grid))
                fac = rbind( fac[2:(2*radius+1),], t(c( rep(NA,radius), fac0, rep(NA,radius) )) )

            # Empty results matrix:
            res = matrix(NA, ncol = in.hdr$ncol, nrow = N.out)

            # Apply the 'fun'ction to each grid column:
            for (j in 1:in.hdr$ncol) {
                wl = as.list(1:N.in)
                for (k in 1:N.in) {
                    wl[[k]] = vl[[k]][,j:(j+2*radius)]
                    if (search.mode == "circle")  wl[[k]][f] = NA
                    if (k == 1) {
                        wk = wl[[1]]
                    } else
                        wk = wk & wl[[k]]
                }
                
                # Use only data from areas within the same zone
                # as defined by the factor grid:
                if (!is.null(in.factor.grid)) {
                    facw = fac[,j:(j+2*radius)]
                    the.fac = centervalue(facw)
                    if (is.na(the.fac)) {
                        wk = wk & NA
                    } else wk[ facw != the.fac ] = NA
                }
                
                if (!all(is.na(wk))) {
                    # Mask NA areas in each of the layers:
                    for (k in 1:N.in) 
                        wl[[k]][ is.na(wk) ] = NA
                        
                    # Convert to vector? Remove NAs?
                    if (mw.to.vector) {
                        for (k in 1:N.in) {
                            wl[[k]] = as.vector(wl[[k]])
                            if (mw.na.rm) 
                                wl[[k]] = wl[[k]][!is.na(wl[[k]])]
                        }
                    }
                    
                    # Pass the (x,y) coordinates to the function?
                    if (pass.location) {
                        x.coord = in.hdr$xllcenter + (j-1) * in.hdr$cellsize
                        loc = list( location = 
                                data.frame(x = x.coord, y = y.coord) )
                    }
                    
                    # Set up list of arguments:
                    names(wl) = in.varnames
                    wl = c( wl, loc, alist(...) )
                    
                    # Call the focal function:
                    r = do.call(fun, wl)

                    # Store the results vector:
                    res[,j] = r
                }
            }
            
            # Replace NA by the no-data value:
            res[ is.na(res) ] = out.nodata.value
            
            # Write one line in each of the output grids:
            for (k in 1:N.out) {
                txt = paste(sprintf(fmt,res[k,]),collapse=" ")
                if (dec!=".") txt = gsub(".",dec,txt,fixed=TRUE)
                writeLines(txt,con = out.files[[k]])
            }
        }
    } # end if (radius > 0)
    
    if (!quiet)  cat("\nDone.\n")

    return(out.filenames)
}




#' Helper function for applying predict methods to stacks of grids.
#'
#' This function can be used to apply the predict method of hopefully any fitted predictive model pixel by pixel to a stack of grids representing the explanatory variables. It is intended to be called primarily by \code{\link{multi.local.function}} or \code{\link{multi.focal.function}}.
#' @name grid.predict
#' @param fit a model object for which prediction is desired
#' @param predfun optional prediction function; if missing, the \code{fit}'s \code{\link{predict}} method is called. In some cases it may be convenient to define a wrapper function for the predict method that may be passed as \code{predfun} argument.
#' @param trafo an optional \code{function(x)} that takes a \code{data.frame} \code{x} and returns a \code{data.frame} with the same number of rows; this is intended to perform transformations on the input variables, e.g. derive a log-transformed variable from the raw input read from the grids, or more complex variables such as the NDVI etc.; the \code{data.frame} resulting from a call to \code{trafo} (if provided) is passed to \code{predfun}
#' @param control.predict an optional list of arguments to be passed on to \code{predfun}; this may be e.g. \code{type="response"} to obtain probability prediction maps from a logistic regression model
#' @param predict.column optional character string: Some predict methods (e.g. \code{predict.lda}) return a data.frame with several columns, e.g. one column per class in a classification problem. \code{predict.column} is used to pick the one that is of interest
#' @param trace integer >=0: positive values give more (=2) or less (=1) information on predictor variables and predictions
#' @param location optional location data received from \code{multi.focal.function}; is added to the \code{newdata} object that is passed on to \code{predfun}.
#' @param \dots these arguments are provided by the calling function, usually \code{\link{multi.local.function}} or \code{\link{multi.focal.function}}.  They contain the explanatory (predictor) variables required by the \code{fit} model.
#' @details \code{grid.predict} is a simple wrapper function. First it binds the arguments in \code{\dots} together in a \code{data.frame} with the raw predictor variables that have been read from their grids by the caller, \code{\link{multi.local.function}} (or \code{\link{multi.focal.function}}). Then it calls the optional \code{trafo} function to transform or combine predictor variables (e.g. perform log transformations, ratioing, arithmetic operations such as calculating the NDVI). Finally the \code{predfun} (or, typically, the default \code{\link{predict}} method of \code{fit}) is called, handing over the \code{fit}, the predictor \code{data.frame}, and the optional \code{control.predict} arguments.
#' @return \code{grid.predict} returns the result of the call to \code{predfun} or the default \code{\link{predict}} method.
#' @references Brenning, A. (2008): Statistical geocomputing combining R and SAGA: The example of landslide susceptibility analysis with generalized additive models. In: J. Boehner, T. Blaschke, L. Montanarella (eds.), SAGA - Seconds Out (= Hamburger Beitraege zur Physischen Geographie und Landschaftsoekologie, 19), 23-32.
#' @author Alexander Brenning
#' @note Though \code{grid.predict} can in principle deal with \code{predict} methods returning factor variables, its usual caller \code{\link{multi.local.function}} / \code{\link{multi.focal.function}} cannot; classification models should be dealt with by setting a \code{type="prob"} (for \code{rpart}) or \code{type="response"} (for logistic regression and logistic additive model) argument, for example (see second Example below).
#' @seealso \code{\link{focal.function}}, \code{\link{multi.local.function}}, \code{\link{multi.focal.function}}
#' @examples 
#' \dontrun{
#' # Assume that d is a data.frame with point observations
#' # of a numerical response variable y and predictor variables
#' # a, b, and c.
#' # Fit a generalized additive model to y,a,b,c.
#' # We want to model b and c as nonlinear terms:
#' require(gam)
#' fit <- gam(y ~ a + s(b) + s(c), data = d)
#' multi.local.function(in.grids = c("a", "b", "c"),
#'     out.varnames = "pred",
#'     fun = grid.predict, fit = fit )
#'     # Note that the 'grid.predict' uses by default the
#'     # predict method of 'fit'.
#' # Model predictions are written to a file named pred.asc
#' }
#' 
#' \dontrun{
#' # A fake example of a logistic additive model:
#' require(gam)
#' fit <- gam(cl ~ a + s(b) + s(c), data = d, family = binomial)
#' multi.local.function(in.grids = c("a", "b", "c"),
#'     out.varnames = "pred",
#'     fun = grid.predict, fit = fit,
#'     control.predict = list(type = "response") )
#'     # 'control.predict' is passed on to 'grid.predict', which
#'     # dumps its contents into the arguments for 'fit''s
#'     # 'predict' method.
#' # Model predictions are written to a file named pred.asc
#' }
#' @keywords spatial
#' @export
grid.predict = function(fit, predfun, trafo, control.predict,
    predict.column, trace = 0, location, ...) 
{
    if (missing(fit)) stop("'fit' object required\n")

    if (trace >= 2 & !missing(location))
        print(utils::str(location))
    
    newdata = as.data.frame( list(...) )

    if (!missing(location)) {
        if (is.vector(location)) # shouldn't be...??
            location = as.data.frame(t(location))
        newdata = cbind(newdata, location)
    }
    
    # Apply transformation function to predictor data.frame:
    if (!missing(trafo))
        newdata = trafo(newdata)

    if (trace >= 2)
        print(utils::str(newdata))
    
    args = list(object = fit, newdata = newdata)
    args = c(args, control.predict)
    
    if (missing(predfun)) {
        pred = do.call( stats::predict, args )
    } else
        pred = do.call( predfun, args )
        
    if (!missing(predict.column))
        pred = pred[,predict.column]
    
    if (trace >= 1)
        print(utils::str(pred))
    
    return(pred)
}




#' @rdname multi.focal.function
#' @name multi.local.function
#' @export
multi.local.function = function( 
    in.grids, in.grid.prefix,
    out.grid.prefix,
    path = NULL, in.path = path, out.path = path,
    fun, in.varnames, out.varnames,
    na.strings = "NA",
    valid.ranges, nodata.values = c(), out.nodata.value, 
    digits = 4, hdr.digits = 10, dec = ".", quiet = TRUE, nlines = Inf,
    na.action = stats::na.exclude, pass.location = FALSE, 
    ... )
{
    # build input filenames:
    if (missing(in.grid.prefix)) in.grid.prefix = ""
    if (is.null(in.grid.prefix)) in.grid.prefix = ""
    if (missing(in.varnames)) {
        in.varnames = in.grids
    } else if (missing(in.grids)) {
        in.grids = in.varnames
    }
    stopifnot(length(in.varnames) == length(unique(in.varnames)))
    stopifnot(length(in.grids) == length(unique(in.grids)))
    stopifnot(length(in.varnames) == length(in.grids))
    do.paste.in = (in.varnames!="") & (in.grid.prefix!="")
    in.filenames = paste( in.grid.prefix, c("","_")[do.paste.in+1], in.grids, sep="" )
    in.filenames = default.file.extension(in.filenames,".asc")
    if (!is.null(in.path)) if (any(in.path != ""))
        in.filenames = file.path(in.path, in.filenames)

    # prepare input files:
    N.in = length(in.filenames)
    in.files = in.hdrs = nodata.vals = as.list(1:N.in)
    for (k in 1:N.in) {
        in.files[[k]] = file(in.filenames[k],open="r")
        in.hdrs[[k]] = read.ascii.grid.header(in.files[[k]],dec=dec)
        nodata.vals[[k]] = unique(c(nodata.values,in.hdrs[[k]]$nodata_value))
        if (k > 1) {
            if ( in.hdrs[[k]]$cellsize != in.hdrs[[1]]$cellsize |
                  in.hdrs[[k]]$ncols != in.hdrs[[1]]$ncols |
                  in.hdrs[[k]]$nrows != in.hdrs[[1]]$nrows )
                stop("incompatible input grids")
        }
    }
    on.exit( for (k in 1:N.in) close(in.files[[k]]), add = TRUE )
    in.hdr = in.hdrs[[1]]
    nlines = max( 1, min( c(nlines, in.hdr$nrows), na.rm = TRUE ) )

    # build output filenames:
    if (missing(out.varnames)) {
        # check if the function will return a vector with variable names
        # when called without arguments:
        out.varnames = try(do.call(fun,list()), silent = TRUE)
        if (missing(out.varnames) || class(out.varnames) == "try-error") {
            if (is.character(fun)) {
                out.varnames = gsub(".","",fun,fixed=TRUE)
            } else if (is.function(fun)) {
                out.varnames = deparse(substitute(fun))
            } else stop("unable to determine 'out.varnames' from 'fun'")
            out.varnames = abbreviate(out.varnames,6)
        }
    }
    if (missing(out.grid.prefix)) out.grid.prefix = ""
    if (is.null(out.grid.prefix)) out.grid.prefix = ""
    stopifnot(length(out.varnames) == length(unique(out.varnames)))
    do.paste = (out.varnames!="") & (out.grid.prefix!="")
    out.filenames = paste( out.grid.prefix, c("","_")[do.paste+1], 
                    out.varnames, sep="" )
    out.filenames = default.file.extension(out.filenames,".asc")
    if (!is.null(out.path)) if (out.path!="")
        out.filenames = file.path(out.path,out.filenames)
    if (any(out.filenames %in% in.filenames))
        stop("one of the output file names is equal to an input file name")

    # prepare output files:
    N.out = length(out.filenames)
    out.files = as.list(1:N.out)
    out.hdr = in.hdr
    if (missing(out.nodata.value)) out.nodata.value = in.hdr$nodata_value
    out.hdr$nodata_value = out.nodata.value
    for (k in 1:N.out) {
        out.files[[k]] = file(out.filenames[k],open="w")
        write.ascii.grid.header(out.files[[k]],out.hdr,dec=dec,hdr.digits=hdr.digits)
    }
    on.exit( for (k in 1:N.out) close(out.files[[k]]), add=TRUE )
    
    if (missing(valid.ranges)) {
        valid.ranges = list()
        for (k in 1:N.in) valid.ranges[[k]] = c(-Inf, Inf)
    }

    fmt = paste("%.",digits,"f",sep="")
    loc = NULL

    # Apply 'fun' as a local function:
    
        # Process one line at a time:
        for (i in 1:nlines) {
            if (!quiet) if ((i %% 10)==0) cat("*")
            if (!quiet) if ((i %% 100)==0) cat("\n")
            
            y.coord = in.hdr$yllcenter + (in.hdr$nrows - i) * in.hdr$cellsize
            
            # Read one line at a time, file by file:
            vl0 = as.list(1:N.in)
            for (k in 1:N.in) {
                vl0[[k]] = scan(in.files[[k]], nlines = 1, quiet = TRUE, 
                    dec = dec, na.strings = na.strings)
                if (length(vl0[[k]]) != in.hdr$ncols) {
                    warning("grid line does not have NCOLS values")
                    vl0[[k]] = c( vl0[[k]], 
                            rep(NA, in.hdr$ncols - length(vl0[[k]])) )
                }
                for (na in nodata.vals[[k]]) 
                    vl0[[k]][ vl0[[k]] == na ] = NA
                vl0[[k]][ vl0[[k]] < valid.ranges[[k]][1] ] = NA
                vl0[[k]][ vl0[[k]] > valid.ranges[[k]][2] ] = NA
            }
                            
            # Pass the (x,y) coordinates to the function?
            if (pass.location) {
                x.coord = c(0:(in.hdr$ncol-1)) * in.hdr$cellsize + in.hdr$xllcenter
                loc = list( location = data.frame(x = x.coord, y = y.coord) )
            }
              
            # Transfer data into argument list:
            args = as.list(1:N.in)
            for (k in 1:N.in) args[[k]] = vl0[[k]]
            ### can this be avoided: args <- vl0 ???
            
            # Remove missing data:               
            args0 = na.action( as.data.frame(args) )
            names(args0) = in.varnames
            unsel = attr(args0, "na.action")
                
            # Any data left?
            if (nrow(args0) > 0) {
                args0 = c( as.list(args0), loc, alist(...) )
                if (is.null(unsel)) {
                    # No NA's:
                    res = matrix(do.call(fun, args0), ncol = in.hdr$ncols, nrow = N.out)
                } else {
                    # NA's were present --> have to re-assign results to non-NA grid cells:
                    res = matrix(NA, ncol = in.hdr$ncols, nrow = N.out)
                    res[,-unsel] = do.call(fun, args0)
                }
            } else res = matrix(NA, ncol = in.hdr$ncols, nrow = N.out)

            # fill in nodata_values:
            res[ is.na(res) ] = out.nodata.value
            for (k in 1:N.out) {
                # convert to character string and write to file:
                txt = paste( sprintf(fmt, res[k,]), collapse = " " )
                if (dec != ".") txt = gsub(".", dec, txt, fixed = TRUE)
                writeLines(txt,con = out.files[[k]])
            }
            rm(res)
        }
        
    if (!quiet)  cat("\nDone.\n")

    return(out.filenames)
}
