#' Define target grid for interpolation
#'
#' Define the resolution and extent of a target grid for interpolation by SAGA modules based on (1) user-provided x/y coordinates, (2) an existing SAGA grid file, or (3) the header data of an ASCII grid. Intended to be used with RSAGA's interpolation functions.
#' @name rsaga.target
#' @param target character: method used for defining the target grid
#' @param user.cellsize Only for \code{target="user.defined"}: raster resolution (in the grid's map units)
#' @param user.x.extent See \code{user.y.extent}
#' @param user.y.extent Only for \code{target="user.defined"}: numeric vectors of length 2: minimum and maximum coordinates of grid cell center points
#' @param target.grid Only for \code{target="target.grid"}: character string giving the name of a SAGA grid file that specifies the extent and resolution of the target grid; this target grid file may be overwritten, depending on the specifics of the SAGA GIS module used.
#' @param header Only for \code{target="header"}: list: ASCII grid header (as returned e.g. by \code{\link{read.ascii.grid.header}}) or defined manually; must at least have components \code{ncols}, \code{nrows}, \code{cellsize}, and either \code{x/yllcorner} or \code{x/yllcenter}.
#' @param env A SAGA geoprocessing environment, see \code{\link{rsaga.env}}.)
#' @note This function is to be used with RSAGA functions \code{\link{rsaga.inverse.distance}}, \code{\link{rsaga.nearest.neighbour}} and \code{\link{rsaga.modified.quadratic.shephard}}. Note that these are currently only compatible with SAGA GIS 2.0.5 and higher.
#' @seealso \code{\link{read.ascii.grid.header}}
#' @examples
#' \dontrun{
#' # IDW interpolation of attribute "z" from the point shapefile
#' # 'points.shp' to a grid with the same extent and resolution
#' # as the (pre-existing) geology grid:
#' rsaga.inverse.distance("points", "dem", field = "z", maxdist = 1000,
#'     target = rsaga.target(target="target.grid",
#'     target.grid = "geology"))
#' }
#' @keywords spatial interface
#' @export
rsaga.target = function(
    target = c("user.defined", "target.grid", "header"),
    user.cellsize = 100,
    user.x.extent, user.y.extent,
    target.grid, header, env = rsaga.env() )
{
    if(env$version == "2.0.4")
        stop("'rsaga.target' currently doesn't support SAGA GIS version 2.0.4\n")
    
    target = match.arg.ext(target, base = 0, numeric = TRUE)
    
    if (target == 2) {
        stopifnot(missing(user.x.extent) & missing(user.y.extent) & missing(target.grid))
        target = 0
        user.cellsize = header$cellsize
        if (!any(names(header) == "xllcenter"))
            header$xllcenter = header$xllcorner + header$cellsize / 2
        if (!any(names(header) == "yllcenter"))
            header$yllcenter = header$yllcorner + header$cellsize / 2
        user.x.extent = c(header$xllcenter, header$xllcenter + header$cellsize * (header$ncols-1))
        user.y.extent = c(header$yllcenter, header$yllcenter + header$cellsize * (header$nrows-1))
    }

    param = list(TARGET = target)
    
    if (target == 0) {
        param = c(param,
            USER_SIZE = user.cellsize,
            USER_XMIN = min(user.x.extent),
            USER_XMAX = max(user.x.extent),
            USER_YMIN = min(user.y.extent),
            USER_YMAX = max(user.y.extent))
    } else if (target == 1) {
        stopifnot(missing(user.x.extent) & missing(user.y.extent))
        target.grid = default.file.extension(target.grid, ".sgrd")
        param = c(param,
            GRID_GRID = target.grid)
    }
    return(param)
}



########     Module io_grid_gdal    ########



#' Import Grid Files to SAGA grid format using GDAL
#' 
#' These functions provide simple interfaces for reading and writing grids from/to ASCII grids and Rd files. Grids are stored in matrices, their headers in lists.
#' @name rsaga.import.gdal
#' @param in.grid file name of a grid in a format supported by GDAL
#' @param out.grid output SAGA grid file name; defaults to \code{in.grid} with the file extension being removed; file extension should not be specified, it defaults to \code{.sgrd}
#' @param env RSAGA geoprocessing environment created by \code{\link{rsaga.env}}
#' @param ... additional arguments to be passed to \code{rsaga.geoprocessor}
#' @details The GDAL Raster Import module of SAGA imports grid data from various file formats using the Geospatial Data Abstraction Library (GDAL) by Frank Warmerdam.
#' GDAL Versions are specific to SAGA versions:
#' \itemize{
#' \item SAGA 2.0.7 - 2.0.8: GDAL v.1.8.0
#' \item SAGA 2.1.0 - 2.1.1: GDAL v.1.10.0
#' \item SAGA 2.1.2 - 2.2.0: GDAL v.1.11.0
#' \item SAGA 2.2.1 - 2.2.3: GDAL v.2.1.0 dev}
#' More information is available at \url{http://www.gdal.org/}.
#' 
#' If \code{in.grid} has more than one band (e.g. RGB GEOTIFF), then output grids with file names of the form \eqn{in.grid{\_}01.sgrd}{in.grid_01.sgrd}, \eqn{in.grid{\_}02.sgrd}{in.grid_02.sgrd} etc. are written, one for each band.
#' 
#' The following raster formats are currently supported. Last updated for SAGA GIS 2.2.3;
#' for a list for a specific SAGA GIS version call \code{rsaga.html.help("io_gdal","GDAL: Import Raster", env = rsaga.env(path="SAGA_Version_to_Test"))}
#' \itemize{
#' \item BAG - Bathymetry Attributed Grid
#' \item ECW - ERDAS Compressed Wavelets (SDK 3.x)
#' \item JP2ECW - ERDAS JPEG2000 (SDK 3.x)
#' \item FITS - Flexible Image Transport System
#' \item GMT - GMT NetCDF Grid Format
#' \item HDF4 - Hierarchical Data Format Release 4
#' \item HDF4Image - HDF4 Dataset
#' \item HDF5 - Hierarchical Data Format Release 5
#' \item HDF5Image - HDF5 Dataset
#' \item KEA - KEA Image Format (.kea)
#' \item MG4Lidar - MrSID Generation 4 / Lidar (.sid)
#' \item MrSID - Multi-resolution Seamless Image Database (MrSID)
#' \item netCDF - Network Common Data Format
#' \item PostgreSQL - PostgreSQL/PostGIS
#' \item VRT - Virtual Raster
#' \item GTiff - GeoTIFF
#' \item NITF - National Imagery Transmission Format
#' \item RPFTOC - Raster Product Format TOC format
#' \item ECRGTOC - ECRG TOC format
#' \item HFA - Erdas Imagine Images (.img)
#' \item SAR_CEOS - CEOS SAR Image
#' \item CEOS - CEOS Image
#' \item JAXAPALSAR - JAXA PALSAR Product Reader (Level 1.1/1.5)
#' \item GFF - Ground-based SAR Applications Testbed File Format (.gff)
#' \item ELAS - ELAS
#' \item AIG - Arc/Info Binary Grid
#' \item AAIGrid - Arc/Info ASCII Grid
#' \item GRASSASCIIGrid - GRASS ASCII Grid
#' \item SDTS - SDTS Raster
#' \item DTED - DTED Elevation Raster
#' \item PNG - Portable Network Graphics
#' \item JPEG - JPEG JFIF
#' \item MEM - In Memory Raster
#' \item JDEM - Japanese DEM (.mem)
#' \item GIF - Graphics Interchange Format (.gif)
#' \item BIGGIF - Graphics Interchange Format (.gif)
#' \item ESAT - Envisat Image Format
#' \item BSB - Maptech BSB Nautical Charts
#' \item XPM - X11 PixMap Format
#' \item BMP - MS Windows Device Independent Bitmap
#' \item DIMAP - SPOT DIMAP
#' \item AirSAR - AirSAR Polarimetric Image
#' \item RS2 - RadarSat 2 XML Product
#' \item SAFE - Sentinel SAFE Product
#' \item PCIDSK - PCIDSK Database File
#' \item PCRaster - PCRaster Raster File
#' \item ILWIS - ILWIS Raster Map
#' \item SGI - SGI Image File Format 1.0
#' \item SRTMHGT - SRTMHGT File Format
#' \item Leveller - Leveller heightfield
#' \item Terragen - Terragen heightfield
#' \item ISIS3 - USGS Astrogeology ISIS cube (Version 3)
#' \item ISIS2 - USGS Astrogeology ISIS cube (Version 2)
#' \item PDS - NASA Planetary Data System
#' \item VICAR - MIPL VICAR file
#' \item TIL - EarthWatch .TIL
#' \item ERS - ERMapper .ers Labelled
#' \item JP2OpenJPEG - JPEG-2000 driver based on OpenJPEG library
#' \item L1B - NOAA Polar Orbiter Level 1b Data Set
#' \item FIT - FIT Image
#' \item GRIB - GRIdded Binary (.grb)
#' \item RMF - Raster Matrix Format
#' \item WCS - OGC Web Coverage Service
#' \item WMS - OGC Web Map Service
#' \item MSGN - EUMETSAT Archive native (.nat)
#' \item RST - Idrisi Raster A.1
#' \item INGR - Intergraph Raster
#' \item GSAG - Golden Software ASCII Grid (.grd)
#' \item GSBG - Golden Software Binary Grid (.grd)
#' \item GS7BG - Golden Software 7 Binary Grid (.grd)
#' \item COSAR - COSAR Annotated Binary Matrix (TerraSAR-X)
#' \item TSX - TerraSAR-X Product
#' \item COASP - DRDC COASP SAR Processor Raster
#' \item R - R Object Data Store
#' \item MAP - OziExplorer .MAP
#' \item PNM - Portable Pixmap Format (netpbm)
#' \item DOQ1 - USGS DOQ (Old Style)
#' \item DOQ2 - USGS DOQ (New Style)
#' \item ENVI - ENVI .hdr Labelled
#' \item EHdr - ESRI .hdr Labelled
#' \item GenBin - Generic Binary (.hdr Labelled)
#' \item PAux - PCI .aux Labelled
#' \item MFF - Vexcel MFF Raster
#' \item MFF2 - Vexcel MFF2 (HKV) Raster
#' \item FujiBAS - Fuji BAS Scanner Image
#' \item GSC - GSC Geogrid
#' \item FAST - EOSAT FAST Format
#' \item BT - VTP .bt (Binary Terrain) 1.3 Format
#' \item LAN - Erdas .LAN/.GIS
#' \item CPG - Convair PolGASP
#' \item IDA - Image Data and Analysis
#' \item NDF - NLAPS Data Format
#' \item EIR - Erdas Imagine Raw
#' \item DIPEx - DIPEx
#' \item LCP - FARSITE v.4 Landscape File (.lcp)
#' \item GTX - NOAA Vertical Datum .GTX
#' \item LOSLAS - NADCON .los/.las Datum Grid Shift
#' \item NTv2 - NTv2 Datum Grid Shift
#' \item CTable2 - CTable2 Datum Grid Shift
#' \item ACE2 - ACE2
#' \item SNODAS - Snow Data Assimilation System
#' \item KRO - KOLOR Raw
#' \item ROI_PAC - ROI_PAC raster
#' \item ISCE - ISCE raster
#' \item ARG - Azavea Raster Grid format
#' \item RIK - Swedish Grid RIK (.rik)
#' \item USGSDEM - USGS Optional ASCII DEM (and CDED)
#' \item GXF - GeoSoft Grid Exchange Format
#' \item NWT_GRD - Northwood Numeric Grid Format .grd/.tab
#' \item NWT_GRC - Northwood Classified Grid Format .grc/.tab
#' \item ADRG - ARC Digitized Raster Graphics
#' \item SRP - Standard Raster Product (ASRP/USRP)
#' \item BLX - Magellan topo (.blx)
#' \item Rasterlite - Rasterlite
#' \item PostGISRaster - PostGIS Raster driver
#' \item SAGA - SAGA GIS Binary Grid (.sdat)
#' \item KMLSUPEROVERLAY - Kml Super Overlay
#' \item XYZ - ASCII Gridded XYZ
#' \item HF2 - HF2/HFZ heightfield raster
#' \item PDF - Geospatial PDF
#' \item OZI - OziExplorer Image File
#' \item CTG - USGS LULC Composite Theme Grid
#' \item E00GRID - Arc/Info Export E00 GRID
#' \item ZMap - ZMap Plus Grid
#' \item NGSGEOID - NOAA NGS Geoid Height Grids
#' \item MBTiles - MBTiles
#' \item IRIS - IRIS data (.PPI, .CAPPi etc)
#' \item PLMOSAIC - Planet Labs Mosaic
#' \item CALS - CALS (Type 1)
#' \item WMTS - OGC Web Map Tile Service
#' \item ESRI Shapefile - ESRI Shapefile
#' \item MapInfo File - MapInfo File
#' \item UK .NTF - UK .NTF
#' \item OGD_SDTS - SDTS
#' \item S57 - IHO S-57 (ENC)
#' \item DGN - Microstation DGN
#' \item OGR_VRT - VRT - Virtual Datasource
#' \item REC EPIInfo .REC
#' \item Memory - Memory
#' \item BNA - Atlas BNA
#' \item CSV - Comma Separated Value (.csv)
#' \item NAS - NAS - ALKIS
#' \item GML - Geography Markup Language
#' \item GPX - GPX
#' \item LIBKML - Keyhole Markup Language (LIBKML)
#' \item KML - Keyhole Markup Language (KML)
#' \item GeoJSON - GeoJSON
#' \item Interlis 1 - Interlis 1
#' \item Interlis 2 - Interlis 2
#' \item OGR_GMT - GMT ASCII Vectors (.gmt)
#' \item GPKG - GeoPackage
#' \item SQLite - SQLite / Spatialite
#' \item ODBC - ODBC
#' \item WAsP - WAsP .map format
#' \item PGeo - ESRI Personal GeoDatabase
#' \item MSSQLSpatial - Microsoft SQL Server Spatial Database
#' \item MySQL - MySQL
#' \item OpenFileGDB - ESRI FileGDB
#' \item XPlane - X-Plane/Flightgear aeronautical data
#' \item DXF - AutoCAD DXF
#' \item Geoconcept - Geoconcept
#' \item GeoRSS - GeoRSS
#' \item GPSTrackMaker - GPSTrackMaker
#' \item VFK - Czech Cadastral Exchange Data Format
#' \item PGDUMP - PostgreSQL SQL dump
#' \item OSM - OpenStreetMap XML and PDF
#' \item GPSBabel - GPSBabel
#' \item SUA - Tim Newport-Peace's Special Use Airspace Format
#' \item OpenAir - OpenAir
#' \item OGR_PDS - Planetary Data Systems TABLE
#' \item WFS - OGC WFS (Web Feature Service)
#' \item HTF - Hydrographic Transfer Vector
#' \item AeronavFAA - Aeronav FAA
#' \item Geomedia - Geomedia .mdb
#' \item EDIGEO - French EDIGEO exchange format
#' \item GFT - Google Fusion Tables
#' \item GME - Google Maps Engine
#' \item SVG - Scalable Vector Graphics
#' \item CouchDB - CouchDB / GeoCouch
#' \item Cloudant - Cloudant / CouchDB
#' \item Idrisi - Idrisi Vector (.vct)
#' \item ARCGEN - Arc/Info Generate
#' \item SEGUKOOA - SEG-P1 / UKOOA P1/90
#' \item SEG-Y - SEG-Y
#' \item ODS - Open Document/ LibreOffice / OpenOffice Spreadsheet
#' \item XLSX - MS Office Open XML spreadsheet
#' \item ElasticSearch - Elastic Search
#' \item Walk - Walk
#' \item CartoDB - CartoDB
#' \item SXF - Storage and eXchange Format
#' \item Selafin - Selafin
#' \item JML - OpenJUMP JML
#' \item PLSCENES - Planet Labs Scenes API
#' \item CSW - OGC CSW (Catalog Search for the Web)
#' \item IDF - INTREST Data Format
#' \item TIGER - U.S. Census TIGER/Line
#' \item AVCBin - Arc/Info Binary Coverage
#' \item AVCE00 - Arc/Info E00 (ASCII) Coverage
#' \item HTTP - HTTP Fetching Wrapper
#' }
#' @references GDAL website: \url{http://www.gdal.org/}
#' @author Alexander Brenning (R interface), Olaf Conrad / Andre Ringeler (SAGA module), Frank Warmerdam (GDAL)
#' @seealso \code{read.ascii.grid}, \code{rsaga.esri.to.sgrd}, \code{read.sgrd}, \code{read.Rd.grid}
#' @keywords spatial interface file
#' @export
rsaga.import.gdal = function( in.grid, out.grid, env = rsaga.env(), ... )
{
    if (missing(out.grid)) {
        out.grid = set.file.extension(in.grid, "")
        out.grid = substr(out.grid, 1, nchar(out.grid) - 1)
    }
    if (env$version == "2.0.4") {
        param = list( GRIDS = out.grid, FILE = in.grid )
    } else {
        param = list( GRIDS = out.grid, FILES = in.grid )
    }
    
    # Module name change with SAGA 2.2.3
    module = "GDAL: Import Raster"
    if (env$version == "2.2.3"){
        module = "Import Raster"
    }
    
    rsaga.geoprocessor("io_gdal", module = module, 
        param = param, env = env, ...)
}




########       Module io_grid       ########


#' Convert ESRI ASCII/binary grids to SAGA grids
#' 
#' \code{rsaga.esri.to.sgrd} converts grid files from ESRI's ASCII (.asc) and binary (.flt) format to SAGA's (version 2) grid format (.sgrd).
#' @name rsaga.esri.to.sgrd
#' @param in.grids character vector of ESRI ASCII/binary grid files (default file extension: \code{.asc}); files should be located in folder \code{in.path}
#' @param out.sgrds character vector of output SAGA grid files; defaults to \code{in.grids} with file extension being replaced by \code{.sgrd}, which is also the default extension if file names without extension are specified; files will be placed in the current SAGA workspace (default: \code{\link{rsaga.env}()$workspace}, or \code{env$workspace} if an \code{env} argument is provided
#' @param in.path folder with \code{in.grids}
#' @param ... optional arguments to be passed to \code{\link{rsaga.geoprocessor}}, including the \code{env} RSAGA geoprocessing environment
#' @return The type of object returned depends on the \code{intern} argument passed to the \code{\link{rsaga.geoprocessor}}. For \code{intern=FALSE} it is a numerical error code (0: success), or otherwise (default) a character vector with the module's console output.
#' 
#' If multiple \code{in.grids} are converted, the result will be a vector of numerical error codes of the same length, or the combination of the console outputs with \code{c()}.
#' @author Alexander Brenning (R interface), Olaf Conrad (SAGA module)
#' @note This function uses module 1 from the SAGA library \code{io_grid}.
#' @seealso \code{\link{rsaga.esri.wrapper}} for an efficient way of applying RSAGA to ESRI ASCII/binary grids; \code{\link{rsaga.env}} 
#' @keywords spatial interface file
#' @export
rsaga.esri.to.sgrd = function( in.grids, 
    out.sgrds=set.file.extension(in.grids,".sgrd"), in.path, ... )
{
    in.grids = default.file.extension(in.grids,".asc")
    out.sgrds = default.file.extension(out.sgrds,".sgrd")
    if (!missing(in.path))
        in.grids = file.path(in.path,in.grids)
    if (length(in.grids) != length(out.sgrds))
        stop("must have the same number of input and outpute grids")
    res = c()
    for (i in 1:length(in.grids))
        res = c(res, rsaga.geoprocessor("io_grid", "Import ESRI Arc/Info Grid",
            list(FILE=in.grids[i],GRID=out.sgrds[i]),...) )
    invisible(res)
}



#' Convert SAGA grids to ESRI ASCII/binary grids
#' 
#' \code{rsaga.sgrd.to.esri} converts grid files from SAGA's (version 2) grid format (.sgrd) to ESRI's ASCII (.asc)  and binary (.flt) format.
#' @name rsaga.sgrd.to.esri
#' @param in.sgrds character vector of SAGA grid files (\code{.sgrd}) to be converted;  files are expected to be found in folder \code{\link{rsaga.env}()$workspace}, or, if an optional \code{env} argument is provided, in \code{env$workspace}
#' @param out.grids character vector of ESRI ASCII/float output file names; defaults to \code{in.sgrds} with the file extension being replaced by \code{.asc} or \code{.flt}, depending on \code{format}. Files will be placed in folder \code{out.path}, existing files will be overwritten 
#' @param out.path folder for \code{out.grids}
#' @param format output file format, either \code{"ascii"} (default; equivalent: \code{format=1}) for ASCII grids or \code{"binary"} (equivalent: \code{0}) for binary ESRI grids (\code{.flt}). 
#' @param georef character: \code{"corner"} (equivalent numeric code: \code{0}) or \code{"center"} (default; equivalent: \code{1}). Determines whether the georeference will be related to the center or corner of its extreme lower left grid cell.
#' @param prec number of digits when writing floating point values to ASCII grid files; either a single number (to be replicated if necessary), or a numeric vector of length \code{length(in.grids)}
#' @param ... optional arguments to be passed to \code{\link{rsaga.geoprocessor}}, including the \code{env} RSAGA geoprocessing environment
#' @return The type of object returned depends on the \code{intern} argument passed to the \code{\link{rsaga.geoprocessor}}. For \code{intern=FALSE} it is a numerical error code (0: success), or otherwise (default) a character vector with the module's console output.
#' @author Alexander Brenning (R interface), Olaf Conrad (SAGA module)
#' @note This function uses module 0 from the SAGA library \code{io_grid}.
#' @seealso \code{\link{rsaga.esri.wrapper}} for an efficient way of applying RSAGA to ESRI ASCII/binary grids; \code{\link{rsaga.env}}
#' @keywords spatial interface file
#' @export
rsaga.sgrd.to.esri = function( in.sgrds, out.grids, out.path,
    format="ascii", georef="corner", prec=5, ... )
{
    in.sgrds = default.file.extension(in.sgrds,".sgrd")
    format = match.arg.ext(format,choices=c("binary","ascii"),base=0,ignore.case=TRUE,numeric=TRUE)
    georef = match.arg.ext(georef,choices=c("corner","center"),base=0,ignore.case=TRUE,numeric=TRUE)
    if (missing(out.grids))
        out.grids = set.file.extension(in.sgrds, c(".flt",".asc")[format+1])
    out.grids = default.file.extension(out.grids, c(".flt",".asc")[format+1])
    if (!missing(out.path))
        out.grids = file.path(out.path,out.grids)
    if (length(out.grids) != length(in.sgrds))
        stop("must have the same number of input and outpute grids")
    if ((length(prec)==1) & (length(in.sgrds)>1))
        prec = rep(prec,length(in.sgrds))
    if (length(prec) != length(in.sgrds))
        stop("must have same number of in-/output grids and 'prec' parameters (or length(prec)==1)")
    res = c()
    for (i in 1:length(in.sgrds))
        res = c(res, rsaga.geoprocessor("io_grid", "Export ESRI Arc/Info Grid",
            list( GRID=in.sgrds[i], FILE=out.grids[i], FORMAT=format, GEOREF=georef, PREC=prec[i]),
            ...))
    invisible(res)
}


#
########    Module ta_morphometry   ########

#' Slope, Aspect, Curvature
#' 
#' Calculates local morphometric terrain attributes (i.e. slope, aspect, and curvatures). Intended for use with SAGA v 2.1.1+. For older versions use \code{\link{rsaga.local.morphometry}}.
#' @name rsaga.slope.asp.curv
#' @param in.dem input: digital elevation model as SAGA grid file (\code{.sgrd})
#' @param out.slope optional output: slope
#' @param out.aspect optional output: aspect
#' @param out.cgene optional output: general curvature
#' @param out.cprof optional output: profile curvature (vertical curvature; degrees)
#' @param out.cplan optional output: plan curvature (horizontal curvature; degrees)
#' @param out.ctang optional output: tangential curvature (degrees)
#' @param out.clong optional output: longitudinal curvature (degrees) Zevenbergen & Thorne (1987) refer to this as profile curvature
#' @param out.ccros optional output: cross-sectional curvature (degrees) Zevenbergen & Thorne (1987) refer to this as the plan curvature
#' @param out.cmini optional output: minimal curvature (degrees)
#' @param out.cmaxi optional output: maximal curvature (degrees)
#' @param out.ctota optional output: total curvature (degrees)
#' @param out.croto optional output: flow line curvature (degrees)
#' @param method character algorithm (see References):
#' \itemize{
#' \item [0] Maximum Slope - Travis et al. (1975) (\code{"maxslope"})
#' \item [1] Max. Triangle Slope - Tarboton (1997) (\code{"maxtriangleslope"})
#' \item [2] Least Squares Fit Plane - Costa-Cabral & Burgess (1996) (\code{"lsqfitplane"})
#' \item [3] Fit 2nd Degree Polynomial - Evans (1979) (\code{"poly2evans"})
#' \item [4] Fit 2nd Degree Polynomial - Heerdegen and Beran (1982) (\code{"poly2heerdegen"})
#' \item [5] Fit 2nd Degree Polynomial - Bauer et al. (1985) (\code{"poly2bauer"})
#' \item [6] default: Fit 2nd Degree Polynomial - Zevenbergen & Thorne (1987) (\code{"poly2zevenbergen"})
#' \item [7] Fit 3rd Degree Polynomial - Haralick (1983) (\code{"poly3haralick"})}
#' @param unit.slope character or numeric (default \code{"radians"}):
#' \itemize{
#' \item [0] \code{"radians"}
#' \item [1] \code{"degrees"}
#' \item [2] \code{"percent"}}
#' @param unit.aspect character or numeric (default is 0, or \code{"radians"}):
#' \itemize{
#' \item [0] \code{"radians"}
#' \item [1] \code{"degrees"}}
#' @param env list, setting up a SAGA geoprocessing environment as created by \code{\link{rsaga.env}}
#' @param ... further arguments to \code{\link{rsaga.geoprocessor}}
#' @details Profile and plan curvature calculation (\code{out.cprof}, \code{out.cplan}) changed in SAGA GIS 2.1.1+ compared to earlier versions. See the following thread on sourceforge.net for an ongoing discussion: \url{http://sourceforge.net/p/saga-gis/discussion/354013/thread/e9d07075/#5727}
#' @return The type of object returned depends on the \code{intern} argument passed to the \code{\link{rsaga.geoprocessor}}. For \code{intern=FALSE} it is a numerical error code (0: success), or otherwise (default) a character vector with the module's console output.
#' @references General references:
#'
#' Jones KH (1998) A comparison of algorithms used to compute hill slope as a property of the DEM. Computers and Geosciences. 24 (4): 315-323.
#' 
#' References on specific methods:
#' 
#' Maximum Slope:
#' 
#' Travis, M.R., Elsner, G.H., Iverson, W.D., Johnson, C.G. (1975): VIEWIT: computation of seen areas, slope, and aspect for land-use planning. USDA F.S. Gen. Tech. Rep. PSW-11/1975, 70 p. Berkeley, California, U.S.A.
#'
#' Maximum Triangle Slope:
#' 
#' Tarboton, D.G. (1997): A new method for the determination of flow directions and upslope areas in grid digital elevation models. Water Ressources Research, 33(2): 309-319.
#'
#' Least Squares or Best Fit Plane:
#'
#' Beasley, D.B., Huggins, L.F. (1982): ANSWERS: User's manual. U.S. EPA-905/9-82-001, Chicago, IL, 54 pp.
#'
#' Costa-Cabral, M., Burges, S.J. (1994): Digital Elevation Model Networks (DEMON): a model of flow over hillslopes for computation of contributing and dispersal areas. Water Resources Research, 30(6): 1681-1692.
#'
#' Fit 2nd Degree Polynomial:
#' 
#' Evans, I.S. (1979): An integrated system of terrain analysis and slope mapping. Final Report on grant DA-ERO-591-73-G0040. University of Durham, England.
#' 
#' Bauer, J., Rohdenburg, H., Bork, H.-R. (1985): Ein Digitales Reliefmodell als Vorraussetzung fuer ein deterministisches  Modell der Wasser- und Stoff-Fluesse. Landschaftsgenese und Landschaftsoekologie, H. 10, Parameteraufbereitung fuer deterministische Gebiets-Wassermodelle, Grundlagenarbeiten zur Analyse von Agrar-Oekosystemen, eds.: Bork, H.-R., Rohdenburg, H., p. 1-15.
#'
#' Heerdegen, R.G., Beran, M.A. (1982): Quantifying source areas through land surface curvature. Journal of Hydrology, 57.
#'
#' Zevenbergen, L.W., Thorne, C.R. (1987): Quantitative analysis of land surface topography. Earth Surface Processes and Landforms, 12: 47-56.
#'
#' Fit 3.Degree Polynomial:
#'
#' Haralick, R.M. (1983): Ridge and valley detection on digital images. Computer Vision, Graphics and Image Processing, 22(1): 28-38.
#'
#' For a discussion on the calculation of slope by ArcGIS check these links:
#'
#' \url{http://forums.esri.com/Thread.asp?c=93&f=1734&t=239914}
#'
#' \url{http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?topicname=how_slope_works}
#' @author Alexander Brenning and Donovan Bangs (R interface), Olaf Conrad (SAGA module)
#' @seealso \code{\link{rsaga.local.morphometry}}, \code{\link{rsaga.parallel.processing}}, \code{\link{rsaga.geoprocessor}},  \code{\link{rsaga.env}}
#' @examples
#' \dontrun{
#' # Simple slope, aspect, and general curvature in degrees:
#' rsaga.slope.asp.curv("lican.sgrd", "slope", "aspect", "curvature",
#'                      method = "maxslope", unit.slope = "degrees", unit.aspect = "degrees")
#' # same for ASCII grids (default extension .asc):
#' rsaga.esri.wrapper(rsaga.slope.asp.curv,
#'                    in.dem="lican", out.slope="slope",
#'                    out.aspect = "aspect", out.cgene = "curvature",
#'                    method="maxslope", unit.slope = "degrees", unit.aspect = "degrees")
#' }
#' @keywords spatial interface
#' @export
rsaga.slope.asp.curv = function(in.dem,
                              out.slope, out.aspect, out.cgene,
                              out.cprof, out.cplan, out.ctang,
                              out.clong, out.ccros, out.cmini,
                              out.cmaxi, out.ctota, out.croto,
                              method = "poly2zevenbergen", 
                              unit.slope = "radians", unit.aspect = "radians",
                              env = rsaga.env(), ...) {
  
  if(env$version != "2.1.1" & env$version != "2.1.2" &
     env$version != "2.1.3" & env$version != "2.1.4" &
     env$version != "2.2.0" & env$version != "2.2.1" &
     env$version != "2.2.2" & env$version != "2.2.3") {
    stop("rsaga.slope.asp.curv only for SAGA GIS 2.1.1+;\n",
         "use rsaga.local.morphometry for older versions of SAGA GIS")
  }
  
  in.dem = default.file.extension(in.dem, ".sgrd")
  method.choices = c("maxslope","maxtriangleslope","lsqfitplane", "poly2evans",
                     "poly2bauer","poly2heerdegen","poly2zevenbergen","poly3haralick")
  if(is.numeric(method) == TRUE)
    stop("Numeric 'method' argument not supported with SAGA GIS 2.1.1+;\n",
         "Use character name of methods - see help(rsaga.slope.asp.curv) for options")
  method = match.arg.ext(method, method.choices, numeric=TRUE, base=0)
  
  unit.slope.choices = c("radians", "degrees", "percent")
  unit.slope = match.arg.ext(unit.slope, unit.slope.choices, numeric=TRUE, base=0)
  
  unit.aspect.choices = c("radians", "degrees")
  unit.aspect = match.arg.ext(unit.aspect, unit.aspect.choices, numeric=TRUE, base=0)
  
  if (missing(out.aspect)) {
    out.aspect = tempfile()
    on.exit(unlink(paste(out.aspect,".*",sep="")), add = TRUE)
  }
  if (missing(out.slope)) {
    out.slope = tempfile()
    on.exit(unlink(paste(out.slope,".*",sep="")), add = TRUE)
  }
  
  param = list(ELEVATION=in.dem, SLOPE=out.slope, ASPECT = out.aspect)
  if(!missing(out.cgene))
    param = c(param, C_GENE = out.cgene)
  if(!missing(out.cprof))
    param = c(param, C_PROF = out.cprof)
  if(!missing(out.cplan))
    param  =c(param, C_PLAN = out.cplan)
  if(!missing(out.ctang))
    param = c(param, C_TANG = out.ctang)
  if(!missing(out.clong))
    param = c(param, C_LONG = out.clong)
  if(!missing(out.ccros))
    param = c(param, C_CROS = out.ccros)
  if(!missing(out.cmini))
    param = c(param, C_MINI = out.cmini)
  if(!missing(out.cmaxi))
    param = c(param, C_MAXI = out.cmaxi)
  if(!missing(out.ctota))
    param = c(param, C_TOTA = out.ctota)
  if(!missing(out.croto))
    param = c(param, C_ROTO = out.croto)
  
  param = c(param, METHOD=method, UNIT_SLOPE=unit.slope, UNIT_ASPECT=unit.aspect)
  
  module = "Slope, Aspect, Curvature"
  
  rsaga.geoprocessor("ta_morphometry", module, param, env = env, ...)
  
  if (!missing(out.cprof) | !missing(out.cplan))
    warning("Plan and profile curvature calculations have changed with SAGA 2.1.1+\n",
            "See help(rsaga.slope.asp.curv) for more information")
}


#' Local Morphometry
#' 
#' Calculates local morphometric terrain attributes (i.e. slope, aspect and curvatures). Intended for use with SAGA versions 2.1.0 and older. Use \code{\link{rsaga.slope.asp.curv}} for SAGA 2.1.1+
#' @name rsaga.local.morphometry
#' @param in.dem input: digital elevation model (DEM) as SAGA grid file (default file extension: \code{.sgrd})
#' @param out.slope optional output: slope (in radians)
#' @param out.aspect optional output: aspect (in radians; north=0, clockwise angles)
#' @param out.curv optional output: curvature
#' @param out.hcurv optional output: horizontal curvature (plan curvature)
#' @param out.vcurv optional output: vertical curvature (profile curvature)
#' @param method character (or numeric): algorithm (see References):
#' \itemize{
#' \item [0] Maximum Slope - Travis et al. (1975) (\code{"maxslope"}, or 0)
#' \item [1] Max. Triangle Slope - Tarboton (1997) (\code{"maxtriangleslope"}, or 1)
#' \item [2] Least Squares Fit Plane - Costa-Cabral and Burgess (1996) (\code{"lsqfitplane"}, or 2)
#' \item [3] Fit 2nd Degree Polynomial - Bauer et al. (1985) (\code{"poly2bauer"}, or 3)
#' \item [4] Fit 2nd Degree Polynomial - Heerdegen and Beran (1982) (\code{"poly2heerdegen"}, or 4)
#' \item [5] default: Fit 2nd Degree Polynomial - Zevenbergen and Thorne (1987) (\code{"poly2zevenbergen"}, or 5)
#' \item [6] Fit 3rd Degree Polynomial - Haralick (1983) (\code{"poly3haralick"}, or 6).}
#' @param env list, setting up a SAGA geoprocessing environment as created by \code{\link{rsaga.env}}
#' @param ... further arguments to \code{\link{rsaga.geoprocessor}}
#' @return The type of object returned depends on the \code{intern} argument passed to the \code{\link{rsaga.geoprocessor}}. For \code{intern=FALSE} it is a numerical error code (0: success), or otherwise (default) a character vector with the module's console output.
#' @references For references and algorithm changes in SAGA GIS 2.1.1+ see \code{\link{rsaga.slope.asp.curv}}.
#' @author Alexander Brenning and Donovan Bangs (R interface), Olaf Conrad (SAGA module)
#' @seealso \code{\link{rsaga.slope.asp.curv}}, \code{\link{rsaga.parallel.processing}}, \code{\link{rsaga.geoprocessor}},  \code{\link{rsaga.env}}
#' @examples
#' \dontrun{
#' # a simple slope algorithm:
#' rsaga.slope("lican.sgrd","slope","maxslope")
#' # same for ASCII grids (default extension .asc):
#' rsaga.esri.wrapper(rsaga.slope,in.dem="lican",out.slope="slope",method="maxslope")
#' }
#' @keywords spatial interface
#' @export
rsaga.local.morphometry = function( in.dem, 
    out.slope, out.aspect, out.curv, out.hcurv, out.vcurv,
    method = "poly2zevenbergen", env = rsaga.env(), ...)
{
  if (!(env$version %in% c("2.0.4","2.0.5","2.0.6","2.0.7","2.0.8","2.0.9","2.1.0"))) {
    rsaga.slope.asp.curv( in.dem=in.dem, out.slope=out.slope, out.aspect=out.aspect, 
        out.cgene=out.curv, out.cplan=out.hcurv, out.cprof=out.vcurv, 
        method=method, env=env, ... )
    warning("rsaga.local.morphometry specific to SAGA versions < 2.1.1\n",
            "Translating provided arguments and using rsaga.slope.asp.curv\n",
            "Note: order of numeric methods have changed with SAGA 2.1.1+")
  } else {
  
    in.dem = default.file.extension(in.dem,".sgrd")
    choices = c("maxslope","maxtriangleslope","lsqfitplane",
        "poly2bauer","poly2heerdegen","poly2zevenbergen","poly3haralick")
    method = match.arg.ext(method,choices,numeric=TRUE,base=0)
    if (missing(out.aspect)) {
        out.aspect = tempfile()
        on.exit(unlink(paste(out.aspect,".*",sep="")), add = TRUE)
    }
    if (missing(out.slope)) {
        out.slope = tempfile()
        on.exit(unlink(paste(out.slope,".*",sep="")), add = TRUE)
    }
    param = list(ELEVATION=in.dem, SLOPE=out.slope, ASPECT=out.aspect)
    if (!missing(out.curv))
        param = c(param, CURV=out.curv)
    if (!missing(out.hcurv))
        param = c(param, HCURV=out.hcurv)
    if (!missing(out.vcurv))
        param = c(param, VCURV=out.vcurv)
    param = c(param, METHOD=method)
    
    module = "Slope, Aspect, Curvature"
    if (any(c("2.0.4","2.0.5","2.0.6") == env$version)) module = "Local Morphometry"
    
    rsaga.geoprocessor("ta_morphometry", module, param, env = env, ...)
  }
    if (!(env$version %in% c("2.0.4","2.0.5","2.0.6","2.0.7","2.0.8","2.0.9","2.1.0"))){
        if (!missing(out.hcurv) | !missing(out.vcurv))
            warning("Plan and profile curvature calculations have changed with SAGA 2.1.1+\n",
                    "See help(rsaga.slope.asp.curv) for more information")
    }
}

#' @rdname rsaga.local.morphometry
#' @name rsaga.slope
#' @export
rsaga.slope = function( in.dem, out.slope, method = "poly2zevenbergen", env = rsaga.env(), ... ) {
    stopifnot(!missing(out.slope))
    if (!(env$version %in% c("2.0.4","2.0.5","2.0.6","2.0.7","2.0.8","2.0.9","2.1.0"))) {
      rsaga.slope.asp.curv( in.dem=in.dem, out.slope=out.slope, method=method, env = env, ... )
    }
    else {
      rsaga.local.morphometry( in.dem=in.dem, out.slope=out.slope, method=method, env = env, ... )
    }
}

#' @rdname rsaga.local.morphometry
#' @name rsaga.aspect
#' @export
rsaga.aspect = function( in.dem, out.aspect, method = "poly2zevenbergen", env = rsaga.env(), ... ) {
    stopifnot(!missing(out.aspect))
    if (!(env$version %in% c("2.0.4","2.0.5","2.0.6","2.0.7","2.0.8","2.0.9","2.1.0"))) {
      rsaga.slope.asp.curv( in.dem=in.dem, out.aspect=out.aspect, method=method, env = env, ... )      
    }
    else {
      rsaga.local.morphometry( in.dem=in.dem, out.aspect=out.aspect, method=method, env = env, ... )
    }
}


#' @rdname rsaga.local.morphometry
#' @name rsaga.curvature
#' @export
rsaga.curvature = function( in.dem, out.curv, method = "poly2zevenbergen", env = rsaga.env(), ... ) {
    stopifnot(!missing(out.curv))
    if (!(env$version %in% c("2.0.4","2.0.5","2.0.6","2.0.7","2.0.8","2.0.9","2.1.0"))) {
      rsaga.slope.asp.curv( in.dem=in.dem, out.cgene=out.curv, method=method, env = env, ... )
    }
    else {
      rsaga.local.morphometry( in.dem=in.dem, out.curv=out.curv, method=method, env = env,  ... )
    }
}

#' @rdname rsaga.local.morphometry
#' @name rsaga.plan.curvature
#' @export
rsaga.plan.curvature = function( in.dem, out.hcurv, method = "poly2zevenbergen", env = rsaga.env(), ... ) {
    stopifnot(!missing(out.hcurv))
    if (!(env$version %in% c("2.0.4","2.0.5","2.0.6","2.0.7","2.0.8","2.0.9","2.1.0"))) {
      rsaga.slope.asp.curv( in.dem=in.dem, out.cplan=out.hcurv, method=method, env = env,  ... )
    }
    else {
      rsaga.local.morphometry( in.dem=in.dem, out.hcurv=out.hcurv, method=method, env = env,  ... )
    }
}

#' @rdname rsaga.local.morphometry
#' @name rsaga.profile.curvature
#' @export
rsaga.profile.curvature = function( in.dem, out.vcurv, method = "poly2zevenbergen", env = rsaga.env(), ... ) {
    stopifnot(!missing(out.vcurv))
    if (!(env$version %in% c("2.0.4","2.0.5","2.0.6","2.0.7","2.0.8","2.0.9","2.1.0"))) {
      rsaga.slope.asp.curv( in.dem=in.dem, out.cprof=out.vcurv, method=method, env = env, ... )
    }
    else {
      rsaga.local.morphometry( in.dem=in.dem, out.vcurv=out.vcurv, method=method, env = env, ... )
    }
}
  

########   Module ta_preprocessor   ########



#' Fill Sinks
#'
#' Several methods for filling closed depressions in digital elevation models that would affect hydrological modeling.
#' @name rsaga.fill.sinks
#' @param in.dem Input: digital elevation model (DEM) as SAGA grid file (default extension: \code{.sgrd}).
#' @param out.dem Output: filled, depression-free DEM (SAGA grid file). Existing files will be overwritten!
#' @param method The depression filling algorithm to be used (character). One of \code{"planchon.darboux.2001"} (default), \code{"wang.liu.2006"}, or \code{"xxl.wang.liu.2006"}.
#' @param out.flowdir (only for \code{"wang.liu.2001"}): Optional output grid file for computed flow directions (see Notes).
#' @param out.wshed (only for \code{"wang.liu.2001"}): Optional output grid file for watershed basins.
#' @param minslope Minimum slope angle (in degree) preserved between adjacent grid cells (default value of \code{0.01} only for \code{method="planchon.darboux.2001"}, otherwise no default).
#' @param ... Optional arguments to be passed to \code{\link{rsaga.geoprocessor}}, including the \code{env} RSAGA geoprocessing environment.
#' @details This function bundles three SAGA modules for filling sinks using three different algorithms (\code{method} argument).
#'
#' \code{"planchon.darboux.2001"}: The algorithm of Planchon and Darboux (2001) consists of increasing the elevation of pixels in closed depressions until the sink disappears and a mininum slope angle of \code{minslope} (default: \code{0.01} degree) is established.
#'
#' \code{"wang.liu.2006"}: This module uses an algorithm proposed by Wang and Liu (2006) to identify and fill surface depressions in DEMs. The method was enhanced to allow the creation of hydrologically sound elevation models, i.e. not only to fill the depressions but also to  preserve a downward slope along the flow path.  If desired, this  is accomplished by preserving a minimum slope gradient (and thus elevation difference) between cells. This is the fully featured version of the module creating a depression-free DEM, a flow path grid and a grid with watershed basins. If you encounter problems processing large data sets (e.g. LIDAR data) with this module try the basic version (\code{xxl.wang.lui.2006}).
#' 
#' \code{"xxl.wang.liu.2006"}: This modified algorithm after Wang and Liu (2006) is designed to work on large data sets.
#' @return The type of object returned depends on the \code{intern} argument passed to the \code{\link{rsaga.geoprocessor}}. For \code{intern=FALSE} it is a numerical error code (0: success), or otherwise (default) a character vector with the module's console output.
#'
#' The function writes SAGA grid files containing of the depression-free preprocessed DEM, and optionally the flow directions and watershed basins.
#' @references Planchon, O., and F. Darboux (2001): A fast, simple and versatile algorithm to fill the depressions of digital elevation models. Catena 46: 159-176.
#'
#' Wang, L. & H. Liu (2006): An efficient method for identifying and filling surface depressions in digital elevation models for hydrologic analysis and modelling. International Journal of Geographical Information Science, Vol. 20, No. 2: 193-213.
#' @author Alexander Brenning (R interface), Volker Wichmann (SAGA module)
#' @note The flow directions are coded as 0 = north, 1 = northeast, 2 = east, ..., 7 = northwest.
#' 
#' If \code{minslope=0}, depressions will only be filled until a horizontal surface is established, which may not be helpful for hydrological modeling.
#' @seealso \code{\link{rsaga.sink.removal}}, \code{\link{rsaga.sink.route}}.
#' @keywords spatial interface
#' @export
rsaga.fill.sinks = function(in.dem,out.dem,
    method="planchon.darboux.2001", out.flowdir, out.wshed, minslope, ...)
{
    stopifnot(is.character(method))
    method = match.arg.ext(method, ignore.case=TRUE, numeric=TRUE, base=2,
        choices=c("planchon.darboux.2001","wang.liu.2006","xxl.wang.liu.2006"))
    in.dem = default.file.extension(in.dem,".sgrd")
    stopifnot(!missing(out.dem))
    if (missing(minslope)) minslope = NULL
    if (method==2) {
        param = list( DEM=in.dem, RESULT=out.dem )
        if (missing(minslope)) minslope = 0.01
        minslope = as.numeric(minslope)
        method = "Fill Sinks (Planchon/Darboux, 2001)"
    } else if (method==3) {
        if (missing(out.flowdir)) {
            out.flowdir = tempfile()
            on.exit(unlink(paste(out.flowdir,".*",sep="")), add = TRUE)
        }
        if (missing(out.wshed)) {
            out.wshed = tempfile()
            on.exit(unlink(paste(out.wshed,".*",sep="")), add = TRUE)
        }
        param = list(ELEV=in.dem, FILLED=out.dem, FDIR=out.flowdir, WSHED=out.wshed)
        method = "Fill Sinks (Wang & Liu)"
    } else if (method==4) {
        param = list(ELEV=in.dem, FILLED=out.dem)
        method = "Fill Sinks XXL (Wang & Liu)"
    }
    if (!is.null(minslope)) param = c( param, MINSLOPE=minslope )
    rsaga.geoprocessor("ta_preprocessor", method, param, ...)
}



#' Sink Drainage Route Detection
#' 
#' Sink drainage route detection.
#' @name rsaga.sink.route
#' @param in.dem input: digital elevation model (DEM) as SAGA grid file (default file extension: \code{.sgrd})
#' @param out.sinkroute output: sink route grid file: non-sinks obtain a value of 0, sinks are assigned an integer between 0 and 8 indicating the direction to which flow from this sink should be routed
#' @param threshold logical: use a threshold value?
#' @param thrsheight numeric: threshold value (default: \code{100})
#' @param ... optional arguments to be passed to \code{\link{rsaga.geoprocessor}}, including the \code{env} RSAGA geoprocessing environment
#' @return The type of object returned depends on the \code{intern} argument passed to the \code{\link{rsaga.geoprocessor}}. For \code{intern=FALSE} it is a numerical error code (0: success), or otherwise (default) a character vector with the module's console output.
#' @author Alexander Brenning (R interface), Olaf Conrad (SAGA module)
#' @note I assume that flow directions are coded as 0 = north, 1 = northeast,  2 = east, ..., 7 = northwest, as in \code{\link{rsaga.fill.sinks}}.
#' @seealso  \code{\link{rsaga.sink.removal}}
#' @examples
#' \dontrun{rsaga.sink.route("dem","sinkroute")
#' rsaga.sink.removal("dem","sinkroute","dem-preproc",method="deepen")}
#' @keywords spatial interface
#' @export
rsaga.sink.route = function(in.dem, out.sinkroute, 
    threshold, thrsheight = 100, ...)
{
    in.dem = default.file.extension(in.dem,".sgrd")
    param = list( ELEVATION=in.dem, SINKROUTE=out.sinkroute )
    if (!missing(threshold)) {
        if (threshold)   param = c( param, THRESHOLD="" )
    }
    # I guess thrsheight is redundant if threshold is missing/false:
    param = c( param, THRSHEIGHT=as.numeric(thrsheight) )
    rsaga.geoprocessor("ta_preprocessor", "Sink Drainage Route Detection", param, ...)
    # was: module = 0
}



#' Sink Removal
#' Remove sinks from a digital elevation model by deepening drainage routes or filling sinks.
#' @name rsaga.sink.removal
#' @param in.dem input: digital elevation model (DEM) as SAGA grid file (default file extension: \code{.sgrd})
#' @param in.sinkroute optional input: sink route grid file
#' @param out.dem output: modified DEM
#' @param method character string or numeric value specifying the algorithm (partial string matching will be applied): \code{"deepen drainage route"} (or 0): reduce the elevation of pixels in order to achieve drainage out of the former sinks \code{"fill sinks"} (or 1): fill sinks until none are left
#' @param ... optional arguments to be passed to \code{\link{rsaga.geoprocessor}}, including the \code{env} RSAGA geoprocessing environment
#' @return The type of object returned depends on the \code{intern} argument passed to the \code{\link{rsaga.geoprocessor}}. For \code{intern=FALSE} it is a numerical error code (0: success), or otherwise (default) a character vector with the module's console output.
#' @author Alexander Brenning (R interface), Olaf Conrad (SAGA module)
#' @note This function uses module 1 from SAGA library \code{ta_preprocessor}.
#' @seealso  \code{\link{rsaga.sink.route}}, \code{\link{rsaga.fill.sinks}}
#' @examples
#' \dontrun{rsaga.sink.route("dem","sinkroute")
#' rsaga.sink.removal("dem","sinkroute","dem-preproc",method="deepen")}
#' @keywords spatial interface
#' @export
rsaga.sink.removal = function(in.dem,in.sinkroute,out.dem,method="fill",...)
{
    in.dem = default.file.extension(in.dem,".sgrd")
    method = match.arg.ext(method,c("deepen drainage routes","fill sinks"),ignore.case=TRUE,numeric=TRUE)
    param = list( DEM=in.dem )
    if (!missing(in.sinkroute)) {
        in.sinkroute = default.file.extension(in.sinkroute,".sgrd")
        param = c(param, SINKROUTE=in.sinkroute)
    }
    param = c( param, DEM_PREPROC=out.dem, METHOD=method )
    rsaga.geoprocessor("ta_preprocessor", "Sink Removal", param, ...)
}





########     Module grid_tools      ########




#' SAGA Modules Close Gaps and Close One Cell Gaps
#'
#' Close (Interpolate) Gaps
#' @name rsaga.close.gaps
#' @param in.dem input: digital elevation model (DEM) as SAGA grid file (default file extension: \code{.sgrd})
#' @param out.dem output: DEM grid file without no-data values (gaps). Existing files will be overwritten!
#' @param threshold tension threshold for adjusting the interpolator (default: 0.1)
#' @param ... optional arguments to be passed to \code{\link{rsaga.geoprocessor}}, including the \code{env} RSAGA geoprocessing environment
#' @details \code{rsaga.close.one.cell.gaps} only fill gaps whose neighbor grid cells have non-missing data.
#'
#' In \code{rsaga.close.gaps}, larger tension thresholds can be used to reduce overshoots and undershoots in the surfaces used to fill (interpolate) the gaps.
#' @return The type of object returned depends on the \code{intern} argument passed to the \code{\link{rsaga.geoprocessor}}. For \code{intern=FALSE} it is a numerical error code (0: success), or otherwise (default) a character vector with the module's console output.
#' @author Alexander Brenning (R interface), Olaf Conrad (SAGA module)
#' @note This function uses modules 7 (\code{rsaga.close.gaps} and 6 \code{rsaga.close.one.cell.gaps} from the SAGA library \code{grid_tools}.
#'
#' SAGA GIS 2.0.5+ has a new additional module \code{Close Gaps with Spline}, which 
#' can be accessed using \code{\link{rsaga.geoprocessor}} (currently no R wrapper 
#' available). See \code{rsaga.get.usage("grid_tools","Close Gaps with Spline")}
#' or in version 2.1.0+ call \code{rsaga.html.help("grid_tools","Close Gaps with Spline")}.
#' @seealso \code{\link{rsaga.geoprocessor}}, \code{\link{rsaga.env}}
#' @examples
#' \dontrun{
#' # using SAGA grids:
#' rsaga.close.gaps("rawdem.sgrd","dem.sgrd")
#' # using ASCII grids:
#' rsaga.esri.wrapper(rsaga.close.gaps,in.dem="rawdem",out.dem="dem")
#' }
#' @keywords spatial interface
#' @export
rsaga.close.gaps = function(in.dem,out.dem,threshold=0.1,...)
{
    in.dem = default.file.extension(in.dem,".sgrd")
    param = list( INPUT=in.dem, RESULT=out.dem, THRESHOLD=as.numeric(threshold) )
    rsaga.geoprocessor("grid_tools", "Close Gaps", param, ...)
}


#' @rdname rsaga.close.gaps
#' @name rsaga.close.one.cell.gaps
#' @keywords spatial interface
#' @export
rsaga.close.one.cell.gaps = function(in.dem,out.dem,...)
{
    in.dem = default.file.extension(in.dem,".sgrd")
    param = list( INPUT = in.dem, RESULT = out.dem )
    rsaga.geoprocessor("grid_tools", "Close One Cell Gaps", 
        param, ...)
}



########     Module ta_lighting     ########



#' Analytical hillshading
#' Analytical hillshading calculation.
#' @name rsaga.hillshade
#' @param in.dem Input digital elevation model (DEM) as SAGA grid file (default extension: \code{.sgrd}).
#' @param out.grid Output hillshading grid (SAGA grid file). Existing files will be overwritten!
#' @param method Available choices (character or numeric): \code{"standard"} (or \code{0} - default), \code{"max90deg.standard"} (\code{1}), \code{"combined.shading"} (\code{2}), \code{"ray.tracing"} (\code{3}). See Details.
#' @param azimuth Direction of the light source, measured in degree  clockwise from the north direction; default 315, i.e. northwest.
#' @param declination Declination of the light source, measured in degree above the horizon (default 45).
#' @param exaggeration Vertical exaggeration of elevation (default: 4). The terrain exaggeration factor allows to increase the shading contrasts in flat areas.
#' @param ... Optional arguments to be passed to \code{\link{rsaga.geoprocessor}}, including the \code{env} RSAGA geoprocessing environment.
#' @details The Analytical Hillshading algorithm is based on the angle between the surface and the incoming light beams, measured in radians.
#' @return The type of object returned depends on the \code{intern} argument passed to the \code{\link{rsaga.geoprocessor}}. For \code{intern=FALSE} it is a numerical error code (0: success), or otherwise (default) a character vector with the module's console output.
#' @author Alexander Brenning (R interface), Olaf Conrad (SAGA module)
#' @note While the default azimuth of 315 degree (northwest) is not physically meaningful on the northern hemisphere, a northwesterly light source is required to properly depict relief in hillshading images. Physically correct southerly light sources results a hillshade that would be considered by most people as inverted: hills look like depressions, mountain chains like troughs.
#' @seealso \code{\link{rsaga.solar.radiation}}, \code{\link{rsaga.insolation}}
#' @examples
#' \dontrun{rsaga.hillshade("dem.sgrd","hillshade")}
#' @keywords spatial interface
#' @export
rsaga.hillshade = function(in.dem, out.grid,
    method="standard", azimuth=315, declination=45, exaggeration=4, ...)
{
    in.dem = default.file.extension(in.dem,".sgrd")
    out.grid = default.file.extension(out.grid,".sgrd")
    method = match.arg.ext(method, numeric=TRUE, ignore.case=TRUE, base=0,
        choices=c("standard","max90deg.standard","combined.shading","ray.tracing"))
    param = list(ELEVATION=in.dem, SHADE=out.grid, METHOD=method,
        AZIMUTH=azimuth, DECLINATION=declination, EXAGGERATION=exaggeration)
    rsaga.geoprocessor("ta_lighting", "Analytical Hillshading", param, ...)
    # was: module = 0
}



#' Potential incoming solar radiation
#'
#' This function calculates the potential incoming solar radiation in an area using different atmospheric models; module available in SAGA GIS 2.0.6+.
#' @name rsaga.pisr
#' @param in.dem name of input digital elevation model (DEM) grid in SAGA grid format (default extension: \code{.sgrd})
#' @param in.svf.grid Optional input grid in SAGA format:  Sky View Factor; see also \code{local.svf} 
#' @param in.vapour.grid Optional input grid in SAGA format:  Water vapour pressure (mbar); see also argument \code{hgt.water.vapour.pressure} 
#' @param in.latitude.grid Optional input grid in SAGA format: Latitude (degree) of each grid cell
#' @param in.longitude.grid see \code{in.latitude.grid}
#' @param out.direct.grid Output grid: Direct insolation (unit selected by \code{unit} argument)
#' @param out.diffuse.grid Output grid: Diffuse insolation
#' @param out.total.grid Optional output grid: Total insolation, i.e. sum of direct and diffuse incoming solar radiation
#' @param out.ratio.grid Optional output grid: Direct to diffuse ratio
#' @param out.duration Optional output grid: Duration of insolation
#' @param out.sunrise Optional output grid: time of sunrise; only calculated if time span is set to single day
#' @param out.sunset Time of sunset; see \code{out.sunrise}
#' @param local.svf logical (default: \code{TRUE}; if TRUE, use sky view factor based on local slope (after Oke, 1988), if no sky view factor grid is provided in \code{in.svf.grid}
#' @param latitude Geographical latitude in degree North (negative values indicate southern hemisphere)
#' @param unit unit of insolation output grids: \code{"kWh/m2"} (default) \code{"kJ/m2"}, or \code{"J/cm2"}
#' @param solconst solar constant, defaults to 1367 W/m2
#' @param enable.bending logical (default: \code{FALSE}): incorporate effects of planetary bending?
#' @param bending.radius Planetary radius, default \code{6366737.96}
#' @param bending.lat.offset if bending is enabled: latitudinal reference  is \code{"user"}-defined (default), or relative to \code{"top"}, \code{"center"} or \code{"bottom"} of grid?
#' @param bending.lat.ref.user user-defined lat. reference for bending, see \code{bending.lat.offset} 
#' @param bending.lon.offset longitudinal reference, i.e. local time,  is \code{"user"}-defined, or relative to \code{"top"}, \code{"center"} (default) or \code{"bottom"} of grid?
#' @param bending.lon.ref.user  user-defined reference for local time (Details??) 
#' @param method specifies how the atmospheric components should be  accounted for: either based on the height of atmosphere and vapour pressure (\code{"height"}, or numeric code 0), or air pressure, water and dust content (\code{"components"}, code 1), or lumped atmospheric transmittance (\code{"lumped"}, code \code{0}) 
#' @param hgt.atmosphere Height of atmosphere (in m); default 12000 m
#' @param hgt.water.vapour.pressure Water vapour pressure in mbar (default 10 mbar); This value is used if no vapour pressure grid is given in  argument \code{in.vapour.grid}
#' @param cmp.pressure atmospheric pressure in mbar, defaults to 1013 mbar
#' @param cmp.water.content water content of a vertical slice of the atmosphere in cm: between 1.5 and 1.7cm, average 1.68cm (default)
#' @param cmp.dust dust factor in ppm; defaults to 100 ppm
#' @param lmp.transmittance transmittance of the atmosphere in percent; usually between 60 (humid areas) and 80 percent (deserts)
#' @param time.range numeric vector of length 2:  time span (hours of the day) for numerical integration
#' @param time.step time step in hours for numerical integration
#' @param start.date list of length two, giving the start date in \code{day} and \code{month} components as numbers; these numbers are one-based (SAGA_CMD uses zero-based numbers internally), i.e. Jan. 1st is \code{list(day=1,month=1)}
#' @param end.date see \code{start.date}
#' @param day.step if \code{days} indicates a range of days, this specifies the time step (number of days) for calculating the incoming solar radiation
#' @param env RSAGA geoprocessing environment obtained with \code{\link{rsaga.env}}; this argument is required for version control (see Note)
#' @param ... optional arguments to be passed to \code{\link{rsaga.geoprocessor}}
#' @details According to SAGA GIS 2.0.7 documentation, "Most options should do well, but TAPES-G based diffuse irradiance calculation ("Atmospheric Effects" methods 2 and 3) needs further revision!" I.e. be careful with \code{method = "components"} and \code{method = "lumped"}.
#' @references 
#' Boehner, J., Antonic, O. (2009): Land surface parameters specific to topo-climatology. In: Hengl, T. and Reuter, H. I. (eds.): Geomorphometry - Concepts, Software, Applications. Elsevier.
#' 
#' Oke, T.R. (1988): Boundary layer climates. London, Taylor and Francis.
#'
#' Wilson, J.P., Gallant, J.C. (eds.), 2000: Terrain analysis - principles and applications. New York, John Wiley and Sons.
#' @author Alexander Brenning (R interface), Olaf Conrad (SAGA module)
#' @note This module is computationally very intensive (depending on the size of the grid and the time resolution, of course). The performance seems to have much improved in SAGA GIS 2.1.0, which by default runs this module in multicore mode (at the release candidate 1 for Windows does).
#'
#' SAGA_CMD uses zero-based days and months, but this R function uses the standard one-based days and months (e.g. day 1 is the first day of the month, month 1 is January) and translates to the SAGA system.
#'
#' This function uses module Potential Incoming Solar Radiation from SAGA library \code{ta_lighting} in SAGA version 2.0.6+.
#' @seealso \code{\link{rsaga.hillshade}}; for similar modules in older SAGA versions (pre-2.0.6) see \code{\link{rsaga.solar.radiation}} and \code{\link{rsaga.insolation}}
#' @keywords spatial interface
#' @export
rsaga.pisr = function(in.dem, in.svf.grid = NULL, in.vapour.grid = NULL, 
    in.latitude.grid = NULL, in.longitude.grid = NULL,
    out.direct.grid, out.diffuse.grid, out.total.grid = NULL, 
    out.ratio.grid = NULL, out.duration, out.sunrise, out.sunset,
    local.svf = TRUE, latitude, 
    unit=c("kWh/m2","kJ/m2","J/cm2"), solconst=1367.0, 
    enable.bending = FALSE, bending.radius = 6366737.96,
    bending.lat.offset = "user", bending.lat.ref.user = 0,
    bending.lon.offset = "center", bending.lon.ref.user = 0,
    method = c("height","components","lumped"),
    hgt.atmosphere = 12000, hgt.water.vapour.pressure = 10,
    cmp.pressure = 1013, cmp.water.content = 1.68, cmp.dust = 100,
    lmp.transmittance = 70,
    time.range = c(0,24), time.step = 0.5,
    start.date = list(day=21, month=3), end.date = NULL, day.step = 5,
    env = rsaga.env(), ...)
{
    if ( (env$version == "2.0.4" | env$version == "2.0.5") ) {
        stop("rsaga.pisr only for SAGA GIS 2.0.6 - 2.2.1;\n",
             " use rsaga.solar.radiation for older versions of SAGA GIS")
    }
    if ( (env$version == "2.2.2" | env$version == "2.2.3") ) {
        stop("rsaga.pisr only for SAGA GIS 2.0.6 - 2.2.1:\n",
             " use rsaga.pisr2 for newer versions of SAGA GIS")
    }

    in.dem = default.file.extension(in.dem,".sgrd")
    if (!is.null(in.svf.grid)) in.svf.grid = default.file.extension(in.svf.grid,".sgrd")
    if (!is.null(in.vapour.grid)) in.vapour.grid = default.file.extension(in.vapour.grid,".sgrd")
    if (!is.null(in.latitude.grid)) in.latitude.grid = default.file.extension(in.latitude.grid,".sgrd")
    if (!is.null(in.longitude.grid)) in.longitude.grid = default.file.extension(in.longitude.grid,".sgrd")
    if (missing(out.direct.grid)) {
        out.direct.grid = tempfile()
        on.exit(unlink(paste(out.direct.grid,".*",sep="")), add = TRUE)
    }
    if (missing(out.diffuse.grid)) {
        out.diffuse.grid = tempfile()
        on.exit(unlink(paste(out.diffuse.grid,".*",sep="")), add = TRUE)
    }
    if (missing(out.total.grid)) {
        out.total.grid = tempfile()
        on.exit(unlink(paste(out.total.grid,".*",sep="")), add = TRUE)
    }
    if (missing(out.ratio.grid)) {
        out.ratio.grid = tempfile()
        on.exit(unlink(paste(out.ratio.grid,".*",sep="")), add = TRUE)
    }
    if (missing(out.duration)) {
        out.duration = tempfile()
        on.exit(unlink(paste(out.duration,".*",sep="")), add = TRUE)
    }
    if (missing(out.sunrise)) {
        out.sunrise = tempfile()
        on.exit(unlink(paste(out.sunrise,".*",sep="")), add = TRUE)
    }
    if (missing(out.sunset)) {
        out.sunset = tempfile()
        on.exit(unlink(paste(out.sunset,".*",sep="")), add = TRUE)
    }

    unit = match.arg.ext(unit,numeric=TRUE,ignore.case=TRUE,base=0)
    method = match.arg.ext(method, numeric = TRUE, ignore.case = TRUE, base = 0)
    bending.lat.offset = match.arg.ext(bending.lat.offset, c("bottom","center","top","user"), 
        numeric = TRUE, ignore.case = TRUE, base = 0)
    bending.lon.offset = match.arg.ext(bending.lon.offset, c("left","center","right","user"), 
        numeric = TRUE, ignore.case = TRUE, base = 0)

    if (!is.null(latitude))
        stopifnot( (latitude>=-90) & (latitude<=90) )
    stopifnot( length(time.range)==2 )
    stopifnot( all(time.range>=0) & all(time.range<=24) & (time.range[1]<time.range[2]) )
    stopifnot( (time.step>0) & (time.step<=12) )
    stopifnot( (day.step>0) & (day.step<=100) )
    stopifnot( is.logical(local.svf) )
    stopifnot( is.logical(enable.bending) )

    param = list( GRD_DEM=in.dem, 
        GRD_DIRECT = out.direct.grid, GRD_DIFFUS = out.diffuse.grid,
        GRD_TOTAL = out.total.grid, GRD_RATIO = out.ratio.grid,
        DURATION = out.duration, 
        SUNRISE = out.sunrise, SUNSET = out.sunset,
        UNITS = unit, SOLARCONST = as.numeric(solconst), LOCALSVF = local.svf,
        BENDING_BENDING = enable.bending,
        METHOD = method,
        #LATITUDE = as.numeric(latitude),  # removed 27 Dec 2011
        DHOUR = time.step )
     
    # Added 27 Dec 2011:   
    if (!is.null(latitude)) {
        stopifnot((latitude >= -90) & (latitude <= 90))
        param = c(param, LATITUDE = as.numeric(latitude))
    }
        
    if (!is.null(in.svf.grid)) param = c( param, GRD_SVF=in.svf.grid )
    if (!is.null(in.vapour.grid)) param = c( param, GRD_VAPOUR=in.vapour.grid )
    stopifnot( !is.null(latitude) | !is.null(in.latitude.grid) ) # added 27 Dec 2011
    if (!is.null(in.latitude.grid)) param = c( param, GRD_LAT=in.latitude.grid )
    if (!is.null(in.longitude.grid)) param = c( param, GRD_LON=in.longitude.grid )

    if (enable.bending) {
        param = c( param,
            BENDING_RADIUS = bending.radius,
            BENDING_LAT_OFFSET = bending.lat.offset,
            BENDING_LAT_REF_USER = bending.lat.ref.user,
            BENDING_LON_OFFSET = bending.lon.offset,
            BENDING_LON_REF_USER = bending.lon.ref.user )
    }
    
    if (method == 0) {
        param = c(param, ATMOSPHERE = as.numeric(hgt.atmosphere),
            VAPOUR = as.numeric(hgt.water.vapour.pressure))
    } else if (method == 1) {
        param = c(param, PRESSURE = as.numeric(cmp.pressure), 
            WATER = as.numeric(cmp.water.content), DUST = as.numeric(cmp.dust))
    } else if (method == 2) {
        stopifnot( (lmp.transmittance>=0) & (lmp.transmittance<=100) )
        param = c(param, LUMPED = as.numeric(lmp.transmittance))
    } else stopifnot( method %in% c(0:2) )
        
    if (is.null(start.date)) { # one year
        stopifnot( is.null(end.date) )
        param = c( param, PERIOD = 2, DAY_A = 0, MONTH_A = 0,
                      DAY_B = 30, MONTH_B = 11 )
    } else {
        if (is.null(end.date)) {
            param = c( param, PERIOD = 1 ) # single day ... or moment (later)
        } else param = c( param, PERIOD = 2 )
        stopifnot(is.list(start.date))
        stopifnot(length(start.date) == 2)
        stopifnot(all(names(start.date %in% c("day","month"))))
        stopifnot( (start.date$day>=1) & (start.date$day<=31) )
        stopifnot( (start.date$month>=1) & (start.date$month<=12) )
        param = c( param, DAY_A = start.date$day - 1,
                    MON_A = start.date$month - 1 )
        if (is.null(end.date)) {
            # check if moment:
            stopifnot(length(time.range) <= 2)
            if (length(time.range) == 2) {
                if (time.range[2] == time.range[1])
                    time.range = time.range[1]
            }
            if (length(time.range) == 1) {
                # moment
                param$PERIOD = 0
                stopifnot(time.range >= 0 & time.range <= 24)
                param = c(param, MOMENT = round(time.range,3))
            } else {
                stopifnot(time.range[1] >= 0 & time.range[1] <= 24)
                stopifnot(time.range[2] >= 0 & time.range[2] <= 24)
                stopifnot(time.range[1] < time.range[2])
                param = c(param, HOUR_RANGE_MIN = time.range[1],
                    HOUR_RANGE_MAX = time.range[2])
            }
        } else {
            # range of days:
            stopifnot(is.list(end.date))
            stopifnot(length(end.date) == 2)
            stopifnot(all(names(end.date %in% c("day","month"))))
            stopifnot( (end.date$day>=1) & (end.date$day<=31) )
            stopifnot( (end.date$month>=1) & (end.date$month<=12) )
            param = c( param, DAY_B = end.date$day - 1,
                        MON_B = end.date$month - 1,
                        DDAYS = day.step )
            if (is.null(time.range)) time.range = c(0,24)
            stopifnot(length(time.range) == 2)
            stopifnot(time.range[1] >= 0 & time.range[1] <= 24)
            stopifnot(time.range[2] >= 0 & time.range[2] <= 24)
            stopifnot(time.range[1] < time.range[2])
            param = c(param, HOUR_RANGE_MIN = time.range[1],
                HOUR_RANGE_MAX = time.range[2])
        }
    }
    
    rsaga.geoprocessor(lib = "ta_lighting", 
        module = "Potential Incoming Solar Radiation",  # = 2
        param = param, env = env, ...)
}

#' Potential incoming solar radiation SAGA 2.2.2+
#'
#' This function calculates the potential incoming solar radiation in an area using different atmospheric models; This function reflects changes to the module with SAGA 2.2.2+.
#' For SAGA versions 2.0.6 to 2.2.1 please see \code{\link{rsaga.pisr}}.
#' @name rsaga.pisr2
#' @param in.dem name of input digital elevation model (DEM) grid in SAGA grid format (default extension: \code{.sgrd})
#' @param in.svf.grid Optional input grid in SAGA format:  Sky View Factor; see also \code{local.svf} 
#' @param in.vapour.grid Optional input grid in SAGA format:  Water vapour pressure (mbar), for use with \code{method = "height"}; default 10 mbar
#' @param in.linke.grid Optional input grid in SAGA format: Linke turbidity coefficient, for use with \code{method = "hofierka"}; default 3.0
#' @param out.direct.grid Output grid: Direct insolation (unit selected by \code{unit} argument)
#' @param out.diffuse.grid Output grid: Diffuse insolation
#' @param out.total.grid Optional output grid: Total insolation, i.e. sum of direct and diffuse incoming solar radiation
#' @param out.ratio.grid Optional output grid: Direct to diffuse ratio
#' @param out.duration Optional output grid: Duration of insolation
#' @param out.sunrise Optional output grid: time of sunrise; only calculated if time span is set to single day
#' @param out.sunset Time of sunset; see \code{out.sunrise}
#' @param local.svf logical (default: \code{TRUE}; if TRUE, use sky view factor based on local slope (after Oke, 1988), if no sky view factor grid is provided in \code{in.svf.grid}
#' @param location specified whether to use constant latitude supplied by \code{latitude} below (\code{"latitude"} or code \code{0}; default) or as calculated from the grid system (\code{"grid"} or code \code{1})
#' @param latitude Geographical latitude in degree North (negative values indicate southern hemisphere)
#' @param unit unit of insolation output grids: \code{"kWh/m2"} (default) \code{"kJ/m2"}, or \code{"J/cm2"}
#' @param solconst solar constant, defaults to 1367 W/m2
#' @param method specifies how the atmospheric components should be  accounted for: either based on the height of atmosphere and vapour pressure (\code{"height"}, or numeric code 0), or air pressure, water and dust content (\code{"components"}, code 1), or lumped atmospheric transmittance (\code{"lumped"}, code \code{2}), or by the method of Hofierka and Suri, 2009 (\code{"hofierka"}, code \code{3}). Default: \code{"lumped"}. 
#' @param hgt.atmosphere Height of atmosphere (in m); default 12000 m. For use with \code{method = "height"}
#' @param cmp.pressure atmospheric pressure in mbar, defaults to 1013 mbar. For use with \code{method = "components"}
#' @param cmp.water.content water content of a vertical slice of the atmosphere in cm: between 1.5 and 1.7cm, average 1.68cm (default). For use with \code{method = "components"}
#' @param cmp.dust dust factor in ppm; defaults to 100 ppm. For use with \code{method = "components"}
#' @param lmp.transmittance transmittance of the atmosphere in percent; usually between 60 (humid areas) and 80 percent (deserts)
#' @param time.range numeric vector of length 2:  time span (hours of the day) for numerical integration
#' @param time.step time step in hours for numerical integration
#' @param start.date list of length three, giving the start date in \code{day}, \code{month}, and \code{year} components as numbers; month is one-based (SAGA_CMD uses zero-based numbers internally), i.e. Jan. 1st 2015 is \code{list(day=1,month=1,year=2015)}
#' @param end.date see \code{start.date}
#' @param day.step if \code{days} indicates a range of days, this specifies the time step (number of days) for calculating the incoming solar radiation
#' @param env RSAGA geoprocessing environment obtained with \code{\link{rsaga.env}}; this argument is required for version control (see Note)
#' @param ... optional arguments to be passed to \code{\link{rsaga.geoprocessor}}
#' @details According to SAGA GIS 2.0.7 documentation, "Most options should do well, but TAPES-G based diffuse irradiance calculation ("Atmospheric Effects" methods 2 and 3) needs further revision!" I.e. be careful with \code{method = "components"} and \code{method = "lumped"}.
#' @references 
#' Boehner, J., Antonic, O. (2009): Land surface parameters specific to topo-climatology. In: Hengl, T. and Reuter, H. I. (eds.): Geomorphometry - Concepts, Software, Applications. Elsevier.
#' 
#' Oke, T.R. (1988): Boundary layer climates. London, Taylor and Francis.
#'
#' Wilson, J.P., Gallant, J.C. (eds.), 2000: Terrain analysis - principles and applications. New York, John Wiley and Sons.
#' 
#' Hofierka, J., Suri, M. (2002): The solar radiation model for Open source GIS: implementation and applications. International GRASS users conference in Trento, Italy, September 2002
#' @author Alexander Brenning & Donovan Bangs (R interface), Olaf Conrad (SAGA module)
#' @note
#' SAGA_CMD uses zero-based months, but this R function uses the standard one-based months (e.g. day 1 is the first day of the month, month 1 is January) and translates to the SAGA system.
#'
#' This function uses module Potential Incoming Solar Radiation from SAGA library \code{ta_lighting} in SAGA version 2.0.6+.
#' Changes to the module with SAGA 2.2.2+ include adding \code{year} to the \code{*.date} arguments to allow calculation across years.
#' The method of Hofierka and Suri (2009) is added, which uses the Linke turbidity coefficient.
#' Duration of insolation (\code{"out.duration"}) is only calculated when the time period is set to a single day.
#' @seealso \code{\link{rsaga.pisr}}; for similar modules in older SAGA versions (pre-2.0.6) see \code{\link{rsaga.solar.radiation}} and \code{\link{rsaga.insolation}}; \code{\link{rsaga.hillshade}}
#' @keywords spatial interface
#' @export
rsaga.pisr2 = function(in.dem, in.svf.grid = NULL, in.vapour.grid = NULL, 
                       in.linke.grid = NULL,
                       out.direct.grid, out.diffuse.grid, out.total.grid = NULL, 
                       out.ratio.grid = NULL, out.duration, out.sunrise, out.sunset,
                       local.svf = TRUE, location = c("latitude", "grid"), latitude = 53, 
                       unit=c("kWh/m2","kJ/m2","J/cm2"), solconst=1367.0,
                       method = c("height","components","lumped","hofierka"),
                       hgt.atmosphere = 12000,
                       cmp.pressure = 1013, cmp.water.content = 1.68, cmp.dust = 100,
                       lmp.transmittance = 70,
                       time.range = c(0,24), time.step = 0.5,
                       start.date = list(day=31, month=10, year=2015), end.date = NULL, day.step = 5,
                       env = rsaga.env(), ...)
{
    if ( env$version != "2.2.2" & env$version != "2.2.3" ) {
        stop("rsaga.pisr2 only for SAGA GIS 2.2.2+;\n",
             " use rsaga.pisr or rsaga.solar.radiation for older versions of SAGA GIS")
    }
    
    in.dem = default.file.extension(in.dem,".sgrd")
    if (!is.null(in.svf.grid)) in.svf.grid = default.file.extension(in.svf.grid,".sgrd")
    if (!is.null(in.vapour.grid)) in.vapour.grid = default.file.extension(in.vapour.grid,".sgrd")
    if (!is.null(in.linke.grid)) in.linke.grid = default.file.extension(in.linke.grid,".sgrd")
    if (missing(out.direct.grid)) {
        out.direct.grid = tempfile()
        on.exit(unlink(paste(out.direct.grid,".*",sep="")), add = TRUE)
    }
    if (missing(out.diffuse.grid)) {
        out.diffuse.grid = tempfile()
        on.exit(unlink(paste(out.diffuse.grid,".*",sep="")), add = TRUE)
    }
    if (missing(out.total.grid)) {
        out.total.grid = tempfile()
        on.exit(unlink(paste(out.total.grid,".*",sep="")), add = TRUE)
    }
    if (missing(out.ratio.grid)) {
        out.ratio.grid = tempfile()
        on.exit(unlink(paste(out.ratio.grid,".*",sep="")), add = TRUE)
    }
    if (missing(out.duration)) {
        out.duration = tempfile()
        on.exit(unlink(paste(out.duration,".*",sep="")), add = TRUE)
    }
    if (missing(out.sunrise)) {
        out.sunrise = tempfile()
        on.exit(unlink(paste(out.sunrise,".*",sep="")), add = TRUE)
    }
    if (missing(out.sunset)) {
        out.sunset = tempfile()
        on.exit(unlink(paste(out.sunset,".*",sep="")), add = TRUE)
    }
    
    unit = match.arg.ext(unit,numeric=TRUE,ignore.case=TRUE,base=0)
    method = match.arg.ext(method, numeric = TRUE, ignore.case = TRUE, base = 0)
    location = match.arg.ext(location, numeric = TRUE, ignore.case = TRUE, base = 0)
    
    if (!is.null(latitude))
        stopifnot( (latitude>=-90) & (latitude<=90) )
    stopifnot( length(time.range)==2 )
    stopifnot( all(time.range>=0) & all(time.range<=24) & (time.range[1]<time.range[2]) )
    stopifnot( (time.step>0) & (time.step<=12) )
    stopifnot( (day.step>0) & (day.step<=100) )
    stopifnot( is.logical(local.svf) )
    
    param = list( GRD_DEM=in.dem, 
                  GRD_DIRECT = out.direct.grid, GRD_DIFFUS = out.diffuse.grid,
                  GRD_TOTAL = out.total.grid, GRD_RATIO = out.ratio.grid,
                  GRD_DURATION = out.duration, 
                  GRD_SUNRISE = out.sunrise, GRD_SUNSET = out.sunset,
                  UNITS = unit, SOLARCONST = as.numeric(solconst), LOCALSVF = local.svf,
                  METHOD = method,
                  HOUR_STEP = time.step )
    
    if (location == 0) {
        if (!is.null(latitude)) {
            stopifnot((latitude >= -90) & (latitude <= 90))
            param = c(param, LATITUDE = as.numeric(latitude))
        }
    } else {
        param = c(param, LOCATION = as.numeric(location))
    }
    
    if (!is.null(in.svf.grid)) param = c( param, GRD_SVF=in.svf.grid )
    if (!is.null(in.vapour.grid)) param = c( param, GRD_VAPOUR=in.vapour.grid )
    if (!is.null(in.linke.grid)) param = c( param, GRD_LINKE=in.linke.grid )
    
    if (method == 0) {
        param = c(param, ATMOSPHERE = as.numeric(hgt.atmosphere))
    } else if (method == 1) {
        param = c(param, PRESSURE = as.numeric(cmp.pressure), 
                  WATER = as.numeric(cmp.water.content), DUST = as.numeric(cmp.dust))
    } else if (method == 2) {
        stopifnot( (lmp.transmittance>=0) & (lmp.transmittance<=100) )
        param = c(param, LUMPED = as.numeric(lmp.transmittance))
    } else if (method == 3) {
        param = param
    } else stopifnot( method %in% c(0:3) )
    
    if (is.null(start.date)) { # one year
        stopifnot( is.null(end.date) )
        param = c( param, PERIOD = 2, DAY_A = 0, MONTH_A = 0,
                   DAY_B = 30, MONTH_B = 11 )
    } else {
        if (is.null(end.date)) {
            param = c( param, PERIOD = 1 ) # single day ... or moment (later)
        } else param = c( param, PERIOD = 2 )
        stopifnot(is.list(start.date))
        stopifnot(length(start.date) == 3)
        stopifnot(all(names(start.date %in% c("day","month","year"))))
        stopifnot( (start.date$day>=1) & (start.date$day<=31) )
        stopifnot( (start.date$month>=1) & (start.date$month<=12) )
        param = c( param, DAY_A = start.date$day ,
                   MON_A = start.date$month - 1,
                   YEAR_A = start.date$year )
        if (is.null(end.date)) {
            # check if moment:
            stopifnot(length(time.range) <= 2)
            if (length(time.range) == 2) {
                if (time.range[2] == time.range[1])
                    time.range = time.range[1]
            }
            if (length(time.range) == 1) {
                # moment
                param$PERIOD = 0
                stopifnot(time.range >= 0 & time.range <= 24)
                param = c(param, MOMENT = round(time.range,3))
            } else {
                stopifnot(time.range[1] >= 0 & time.range[1] <= 24)
                stopifnot(time.range[2] >= 0 & time.range[2] <= 24)
                stopifnot(time.range[1] < time.range[2])
                param = c(param, HOUR_RANGE_MIN = time.range[1],
                          HOUR_RANGE_MAX = time.range[2])
            }
        } else {
            # range of days:
            stopifnot(is.list(end.date))
            stopifnot(length(end.date) == 3)
            stopifnot(all(names(end.date %in% c("day","month","year"))))
            stopifnot( (end.date$day>=1) & (end.date$day<=31) )
            stopifnot( (end.date$month>=1) & (end.date$month<=12) )
            param = c( param, DAY_B = end.date$day,
                       MON_B = end.date$month - 1,
                       YEAR_B = end.date$year,
                       DAYS_STEP = day.step )
            if (is.null(time.range)) time.range = c(0,24)
            stopifnot(length(time.range) == 2)
            stopifnot(time.range[1] >= 0 & time.range[1] <= 24)
            stopifnot(time.range[2] >= 0 & time.range[2] <= 24)
            stopifnot(time.range[1] < time.range[2])
            param = c(param, HOUR_RANGE_MIN = time.range[1],
                      HOUR_RANGE_MAX = time.range[2])
        }
    }
    
    rsaga.geoprocessor(lib = "ta_lighting", 
                       module = "Potential Incoming Solar Radiation",  # = 2
                       param = param, env = env, ...)
}


#' Potential incoming solar radiation
#'
#' This function calculates the potential incoming solar radiation in an area either using a lumped atmospheric transmittance model or estimating it based on water and dust content. Use \code{\link{rsaga.pisr}} instead with SAGA GIS 2.0.6+.
#' @name rsaga.solar.radiation
#' @param in.dem name of input digital elevation model (DEM) grid in SAGA grid format (default extension: \code{.sgrd})
#' @param out.grid output grid file for potential incoming solar radiation sums
#' @param out.duration Optional output grid file for duration of insolation
#' @param latitude Geographical latitude in degree North (negative values indicate southern hemisphere)
#' @param unit unit of the \code{out.grid} output: \code{"kWh/m2"} (default) or \code{"J/m2"}
#' @param solconst solar constant, defaults to 1367 W/m2
#' @param method specifies how the atmospheric components should be accounted for: either based on a lumped atmospheric transmittance as specified by argument \code{transmittance} (\code{"lumped"}, or numeric code \code{0}; default); or by calculating the components corresponding to water and dust (\code{"components"}, code \code{1})
#' @param transmittance transmittance of the atmosphere in percent; usually between 60 (humid areas) and 80 percent (deserts)
#' @param pressure atmospheric pressure in mbar
#' @param water.content water content of a vertical slice of the atmosphere in cm: between 1.5 and 1.7cm, average 1.68cm (default)
#' @param dust dust factor in ppm; defaults to 100ppm
#' @param time.range numeric vector of length 2:  time span (hours of the day) for numerical integration
#' @param time.step time step in hours for numerical integration
#' @param days either a list with components \code{day} and \code{month} specifying a single day of the year for radiation modeling; OR a numeric vector of length 2 specifying the start and end date (see Note below)
#' @param day.step if \code{days} indicates a range of days, this specifies the time step (number of days) for calculating the incoming solar radiation
#' @param env RSAGA geoprocessing environment obtained with \code{\link{rsaga.env}}; this argument is required for version control (see Note)
#' @param ... optional arguments to be passed to \code{\link{rsaga.geoprocessor}}
#' @references Wilson, J.P., Gallant, J.C. (eds.), 2000: Terrain analysis - principles and applications. New York, John Wiley & Sons.
#' @author Alexander Brenning (R interface), Olaf Conrad (SAGA module)
#' @note This module ceased to exist under SAGA GIS 2.0.6+, which has a similar (but more flexible) module Potential Solar Radiation that is interfaced by \code{\link{rsaga.pisr}}.
#'
#' SAGA_CMD uses zero-based days and months, but this R function uses the standard one-based days and months (e.g. day 1 is the first day of the month, month 1 is January) and translates to the SAGA system.
#'
#' In SAGA 2.0.2, solar radiation sums calculated for a range of days, say \code{days=c(a,b)} actually calculate radiation only for days \code{a,...,b-1} (in steps of \code{day.step} - I used \code{day.step=1} in this example).  The setting \code{a=b} however gives the same result as \code{b=a+1}, and indeed \code{b=a+2} gives twice the radiation sums and potential sunshine duration that \code{a=b} and \code{b=a+1} both give.
#'
#' The solar radiation module of SAGA 2.0.1 had a bug that made it impossible to pass a range of \code{days} of the year or a range of hours of the day (\code{time.range}) to SAGA. These options work in SAGA 2.0.1.
#'
#' This function uses module Incoming Solar Radiation from SAGA GIS library \code{ta_lighting}.
#' @seealso \code{\link{rsaga.hillshade}}, \code{\link{rsaga.insolation}}
#' @examples
#' \dontrun{
#' # potential solar radiation on Nov 7 in Southern Ontario...
#' rsaga.solar.radiation("dem","solrad","soldur",latitude=43,
#'     days=list(day=7,month=11),time.step=0.5)
#' }
#' @keywords spatial interface
#' @export
rsaga.solar.radiation = function(in.dem, out.grid, out.duration, latitude, 
    unit=c("kWh/m2","J/m2"), solconst=1367.0, method=c("lumped","components"),
    transmittance=70, pressure=1013, water.content=1.68, dust=100,
    time.range=c(0,24), time.step=1,
    days=list(day=21,month=3), day.step=5,
    env = rsaga.env(), ...)
{
    if ( !(env$version == "2.0.4" | env$version == "2.0.5") ) {
        stop("rsaga.solar.radiation only for SAGA GIS 2.0.4 / 2.0.5;\n",
             " use rsaga.pisr for SAGA GIS 2.0.6+")
    }

    in.dem = default.file.extension(in.dem,".sgrd")
    if (missing(out.duration)) {
        out.duration = tempfile()
        on.exit(unlink(paste(out.duration,".*",sep="")), add = TRUE)
    }
    unit = match.arg.ext(unit,numeric=TRUE,ignore.case=TRUE,base=0)
    method = match.arg.ext(method,numeric=TRUE,ignore.case=TRUE,base=0)
    stopifnot( (transmittance>=0) & (transmittance<=100) )
    stopifnot( (latitude>=-90) & (latitude<=90) )
    stopifnot( length(time.range)==2 )
    stopifnot( all(time.range>=0) & all(time.range<=24) & (time.range[1]<time.range[2]) )
    stopifnot( (time.step>0) & (time.step<=12) )
    stopifnot( (day.step>0) & (day.step<=100) )

    param = list( ELEVATION=in.dem, INSOLAT=out.grid, DURATION=out.duration,
        UNIT=unit, SOLCONST=as.numeric(solconst), METHOD=method,
        TRANSMITT=as.numeric(transmittance), PRESSURE=as.numeric(pressure), 
        WATER=as.numeric(water.content), DUST=as.numeric(dust),
        LATITUDE=as.numeric(latitude), 
        HOUR_RANGE_MIN=time.range[1], HOUR_RANGE_MAX=time.range[2], 
        HOUR_STEP=time.step )
        
    if (is.null(days)) { # one year
        param = c( param, TIMESPAN=2 )
    } else if (is.list(days)) { # single day
        stopifnot(length(days)==2)
        stopifnot( (days$day>=1) & (days$day<=31) )
        stopifnot( (days$month>=1) & (days$month<=12) )
        param = c( param, TIMESPAN=0,
            SINGLE_DAY_DAY=days$day-1, SINGLE_DAY_MONTH=days$month-1 )
    } else if (is.numeric(days)) { # range of days
        stopifnot(length(days)==2)
        stopifnot( days[1] <= days[2] )
        stopifnot( (days[1]>=1) & (days[2]<=366) )
        param = c( param, TIMESPAN=1, 
            DAY_RANGE_MIN=days[1], DAY_RANGE_MAX=days[2], 
            DAY_STEP=day.step )
    }
    rsaga.geoprocessor(lib = "ta_lighting", 
        module = "Incoming Solar Radiation",  # = 2
        param = param, env = env, ...)
}



#' Incoming Solar Radiation (Insolation)
#'
#' This function calculates the amount of incoming solar radiation (insolation) depending on slope, aspect, and atmospheric properties. Module not available in SAGA GIS 2.0.6 and 2.0.7.
#' @name rsaga.insolation
#' @param in.dem Name of input digital elevation model (DEM) grid in SAGA grid format (default extension: \code{.sgrd})
#' @param in.vapour Optional input: SAGA grid file giving the water vapour pressure in mbar
#' @param in.latitude Optional input: SAGA grid file giving for each pixel the latitude in degree
#' @param in.longitude Optional input: SAGA grid file giving for each pixel the longitude in degree
#' @param out.direct Optional output grid file for direct insolation
#' @param out.diffuse Optional output grid file for diffuse insolation
#' @param out.total Optional output grid file for total insolation, i.e. the sum of direct and diffuse insolation
#' @param horizontal logical; project radiation onto a horizontal surface? (default: \code{FALSE}, i.e. use the actual inclined surface as a reference area) 
#' @param solconst solar constant in Joule; default: 8.164 J/cm2/min (=1360.7 kWh/m2; the more commonly used solar constant of 1367 kWh/m2 corresponds to 8.202 J/cm2/min)
#' @param atmosphere height of atmosphere in m; default: 12000m
#' @param water.vapour.pressure if no water vapour grid is given, this argument specifies a constant water vapour pressure that is uniform in space; in mbar, default 10 mbar
#' @param type type of time period: \code{"moment"} (equivalent: \code{0}) for a single instant, \code{"day"} (or \code{1}) for a single day, \code{"range.of.days"} (or \code{2}), or \code{"same.moment.range.of.days"} (or \code{3}) for the same moment in a range of days; default: \code{"moment"} 
#' @param time.step time resolution in hours for discretization within a day
#' @param day.step time resolution in days for a range of days
#' @param days numeric vector of length 2, specifying the first and last day of a range of days (for \code{type}s 2 and 3)
#' @param moment if \code{type="moment"} or \code{"same.moment.range.of.days"}, \code{moment} specifies the time of the day (hour between 0 and 24) for which the insolation is to be calculated
#' @param latitude if no \code{in.latitude} grid is given, this will specify a fixed geographical latitude for the entire grid
#' @param bending should planetary bending be modeled? (default: \code{FALSE})
#' @param radius planetary radius
#' @param lat.offset \code{latitude} relates to grids \code{"bottom"}(equivalent code: \code{0}), \code{"center"} (1), \code{"top"} (2), or \code{"user"}-defined reference (default: \code{"user"}); in the latter case, \code{lat.ref.user} defines the reference
#' @param lat.ref.user if \code{in.latitude} is missing and \code{lat.offset="user"}, then this numeric value defines the latitudinal reference (details??)
#' @param lon.offset local time refers to grid's \code{"left"} edge (code 0), \code{"center"} (1), \code{"right"} edge (2), or a  \code{"user"}-defined reference.
#' @param lon.ref.user if \code{in.longitude} is missing and \code{lon.offset="user"}, then this numeric value defines the reference of the local time (details??)
#' @param ... optional arguments to be passed to \code{\link{rsaga.geoprocessor}}, including the \code{env} RSAGA geoprocessing environment
#' @details Calculation of incoming solar radiation (insolation). Based on the SADO (System for the Analysis of Discrete Surfaces) routines developed  by Boehner & Trachinow.
#' @return The type of object returned depends on the \code{intern} argument passed to the \code{\link{rsaga.geoprocessor}}. For \code{intern=FALSE} it is a numerical error code (0: success), or otherwise (default) a character vector with the module's console output.
#' @author Alexander Brenning (R interface), Olaf Conrad (SAGA module)
#' @note This function uses module \code{Insolation} (code: 3) from SAGA library \code{ta_lighting}. It is availble in SAGA GIS 2.0.4 and 2.0.5 but not 2.0.6 and 2.0.7; see \code{\link{rsaga.pisr}}.
#' @seealso \code{\link{rsaga.solar.radiation}}, \code{\link{rsaga.pisr}},  \code{\link{rsaga.hillshade}}
#' @keywords spatial interface
#' @export
rsaga.insolation = function(in.dem, in.vapour, in.latitude, in.longitude,
    out.direct, out.diffuse, out.total,
    horizontal=FALSE, solconst=8.1640, atmosphere=12000, water.vapour.pressure=10.0,
    type=c("moment","day","range.of.days","same.moment.range.of.days"),
    time.step=1, day.step=5, days, moment, latitude, bending=FALSE,
    radius=6366737.96,
    lat.offset="user", lat.ref.user=0,
    lon.offset="center", lon.ref.user=0,
     ...)
{
    in.dem = default.file.extension(in.dem,".sgrd")
    param = list( GRD_DEM=in.dem )
    type = match.arg.ext(type,numeric=TRUE,ignore.case=TRUE,base=0)
    stopifnot( (!missing(out.direct)) | (!missing(out.diffuse)) | (!missing(out.total)) )
    stopifnot( !missing(latitude) )
    if (!missing(moment)) {
        if (!(type==0 | type==3)) {
            warning("'moment' argument only relevant for 'type=\"moment\"'\n",
                    "or 'type=\"same.moment.range.of.days\"' -\n",
                    "ignoring the 'moment' argument")
        }
    }
    if (!missing(in.vapour)) {
        in.vapour = default.file.extension(in.vapour,".sgrd")
        param = c(param, GRD_VAPOUR=in.vapour)
    }
    if (!missing(in.latitude)) {
        in.latitude = default.file.extension(in.latitude,".sgrd")
        param = c(param, GRD_LAT=in.latitude)
    }
    if (!missing(in.longitude)) {
        in.longitude = default.file.extension(in.longitude,".sgrd")
        param = c(param, GRD_LON=in.longitude)
    }
    if (!missing(out.direct)) param = c(param, GRD_DIRECT=out.direct)
    if (!missing(out.diffuse)) param = c(param, GRD_DIFFUS=out.diffuse)
    if (!missing(out.total)) param = c(param, GRD_TOTAL=out.total)
    stopifnot( (days[1]>=0) & (days[1]<=366) )
    param = c(param, BHORIZON=horizontal, SOLARCONST=solconst,
        ATMOSPHERE=atmosphere, VAPOUR=water.vapour.pressure,
        PERIOD=type, DHOUR=time.step, DDAYS=day.step,
        DAY_A=days[1])
    if (type>=2) { # range of days / same moment in a range of days
        stopifnot( (days[2]>=days[1]) & (days[2]<=366) )
        param = c(param, DAY_B=days[2])
    }
    if ((type==0) | (type==3)) {
        stopifnot( (moment>=0) & (moment<=24) )
        param = c(param, MOMENT=moment)
    }
    param = c(param, LATITUDE=latitude, BENDING=bending, RADIUS=radius)
    lat.offset = match.arg.ext(lat.offset, c("bottom","center","top","user"),
        numeric=TRUE, ignore.case=TRUE, base=0)
    lon.offset = match.arg.ext(lon.offset, c("left","center","right","user"),
        numeric=TRUE, ignore.case=TRUE, base=0)
    param = c(param, LAT_OFFSET=lat.offset)
    if (lat.offset==3) { # user-defined
        #stopifnot(!missing(lat.ref.user))
        param = c(param, LAT_REF_USER=as.numeric(lat.ref.user))
    }
    param = c(param, LON_OFFSET=lon.offset)
    if (lon.offset==3) { # user-defined
        #stopifnot(!missing(lon.ref.user))
        param = c(param, LON_REF_USER=as.numeric(lon.ref.user))
    }
    rsaga.geoprocessor(lib = "ta_lighting", 
        module = "Insolation", # = 3 
        param = param, ...)
}





########     Module grid_filter     ########



#' Simple Filters
#'
#' Apply a smoothing, sharpening or edge filter to a SAGA grid.
#' @name rsaga.filter.simple
#' @param in.grid input: SAGA grid file (default file extension: \code{.sgrd})
#' @param out.grid output: SAGA grid file
#' @param mode character or numeric: shape of moving window, either \code{"square"} (=0) or \code{"circle"} (=1, default)
#' @param method character or numeric: \code{"smooth"} (=0), \code{"sharpen"} (=1), or \code{"edge"} (=2)
#' @param radius positive integer: radius of moving window
#' @param ... optional arguments to be passed to \code{\link{rsaga.geoprocessor}}, including the \code{env} RSAGA geoprocessing environment
#' @return The type of object returned depends on the \code{intern} argument passed to the \code{\link{rsaga.geoprocessor}}. For \code{intern=FALSE} it is a numerical error code (0: success), or otherwise (the default) a character vector with the module's console output.
#' @author Alexander Brenning (R interface), Olaf Conrad (SAGA module)
#' @seealso \code{\link{rsaga.filter.gauss}}
#' @examples \dontrun{rsaga.filter.simple("dem","dem-smooth",radius=4)}
#' @keywords spatial interface
#' @export
rsaga.filter.simple = function(in.grid, out.grid, mode="circle",
    method=c("smooth","sharpen","edge"), radius,...)
{
    in.grid = default.file.extension(in.grid,".sgrd")
    mode = match.arg.ext(mode,choices=c("square","circle"),
        numeric=TRUE,base=0,ignore.case=TRUE)
    method = match.arg.ext(method,numeric=TRUE,base=0,ignore.case=TRUE)
    if (missing(radius)) stop("the search 'radius' argument (in # pixels) must be specified")
    if (round(radius) != radius) {
        warning("'radius' must be an integer >=1 (# pixels); rounding it...")
        radius = round(radius)
    }
    if (radius<1) {
        warning("'radius' must be an integer >=1 (# pixels); setting 'radius=1'...")
        radius = 1
    }
    param = list(INPUT=in.grid, RESULT=out.grid, MODE=mode,
        METHOD=method, RADIUS=radius)
    rsaga.geoprocessor(lib = "grid_filter", 
        module = "Simple Filter",
        param = param, ...)
}



#' Gauss Filter
#' 
#' Smooth a grid using a Gauss filter.
#' @name rsaga.filter.gauss
#' @param in.grid input: SAGA GIS grid file (default file extension: \code{.sgrd})
#' @param out.grid output: SAGA GIS grid file
#' @param sigma numeric, >0.0001: standard deviation parameter of Gauss filter
#' @param radius positive integer: radius of moving window
#' @param ... optional arguments to be passed to \code{\link{rsaga.geoprocessor}}, including the \code{env} RSAGA geoprocessing environment
#' @return The type of object returned depends on the \code{intern} argument passed to the \code{\link{rsaga.geoprocessor}}. For \code{intern=FALSE} it is a numerical error code (0: success), or otherwise (the default) a character vector with the module's console output.
#' @author Alexander Brenning (R interface), Olaf Conrad (SAGA module)
#' @seealso \code{\link{rsaga.filter.simple}}
#' @keywords spatial interface
#' @export
rsaga.filter.gauss = function(in.grid, out.grid, sigma,
    radius=ceiling(2*sigma),...)
{
    in.grid = default.file.extension(in.grid,".sgrd")
    if (missing(sigma)) stop("the 'sigma' standard deviation argument (in # pixels) must be specified")
    stopifnot(sigma>0.0001)
    if (round(radius) != radius) stop("'radius' must be an integer (# pixels)")
    stopifnot(radius>=1)
    param = list(INPUT=in.grid, RESULT=out.grid, SIGMA=sigma, RADIUS=radius)
    rsaga.geoprocessor(lib = "grid_filter", 
        module = "Gaussian Filter", # = 1, 
        param, ...)
}





########     Module ta_hydrology    ########



#' Parallel Processing
#'
#' Calculate the size of the local catchment area (contributing area), the catchment height, catchment slope and aspect, and flow path length, using parallel processing algorithms including the recommended multiple flow direction algorithm. This set of algorithms processes a digital elevation model (DEM) downwards from the highest to the lowest cell.\cr No longer supported with SAGA GIS 2.1.3+. See \code{\link{rsaga.topdown.processing}}.
#' @name rsaga.parallel.processing
#' @param in.dem input: digital elevation model (DEM) as SAGA grid file (default file extension: \code{.sgrd})
#' @param in.sinkroute optional input: SAGA grid with sink routes
#' @param in.weight optional intput: SAGA grid with weights
#' @param out.carea output: catchment area grid
#' @param out.cheight optional output: catchment height grid
#' @param out.cslope optional output: catchment slope grid
#' @param out.caspect optional output: catchment aspect grid
#' @param out.flowpath optional output: flow path length grid
#' @param step integer >=1: step parameter
#' @param method character or numeric: choice of processing algorithm: Deterministic 8 (\code{"d8"} or 0), Rho 8 (\code{"rho8"} or 1), Braunschweiger Reliefmodell (\code{"braunschweig"} or 2), Deterministic Infinity (\code{"dinf"} or 3), Multiple Flow Direction (\code{"mfd"} or 4, the default), Multiple Triangular Flow Direction (\code{"mtfd"}, or 5).
#' @param linear.threshold numeric (number of grid cells): threshold above which linear flow (i.e. the Deterministic 8 algorithm) will be used; linear flow is disabled for \code{linear.threshold=Inf} (the default)
#' @param convergence numeric >=0: a parameter for tuning convergent/ divergent flow; default value of \code{1.1} gives realistic results and should not be changed
#' @param env list, setting up a SAGA geoprocessing environment as created by \code{\link{rsaga.env}}
#' @param ... further arguments to \code{\link{rsaga.geoprocessor}}
#' @details Refer to the references for details on the available algorithms.
#' @return The type of object returned depends on the \code{intern} argument passed to the \code{\link{rsaga.geoprocessor}}. For \code{intern=FALSE} it is a numerical error code (0: success), or otherwise (the default) a character vector with the module's console output.
#' @references
#' Deterministic 8:
#'
#' O'Callaghan, J.F., Mark, D.M. (1984): The extraction of drainage networks from digital elevation data. Computer Vision, Graphics and Image Processing, 28: 323-344.
#'
#' Rho 8:
#'
#' Fairfield, J., Leymarie, P. (1991): Drainage networks from grid digital elevation models. Water Resources Research, 27: 709-717.
#'
#' Braunschweiger Reliefmodell:
#'
#' Bauer, J., Rohdenburg, H., Bork, H.-R. (1985): Ein Digitales Reliefmodell als Vorraussetzung fuer ein deterministisches Modell der Wasser- und Stoff-Fluesse. Landschaftsgenese und Landschaftsoekologie, H. 10, Parameteraufbereitung fuer deterministische Gebiets-Wassermodelle, Grundlagenarbeiten zu Analyse von Agrar-Oekosystemen, eds.: Bork, H.-R., Rohdenburg, H., p. 1-15.
#'
#' Deterministic Infinity:
#'
#' Tarboton, D.G. (1997): A new method for the determination of flow directions and upslope areas in grid digital elevation models. Water Ressources Research, 33(2): 309-319.
#'
#' Multiple Flow Direction:
#'
#' Freeman, G.T. (1991): Calculating catchment area with divergent flow based on a regular grid. Computers and Geosciences, 17: 413-22.
#'
#' Quinn, P.F., Beven, K.J., Chevallier, P., Planchon, O. (1991): The prediction of hillslope flow paths for distributed hydrological modelling using digital terrain models. Hydrological Processes, 5: 59-79.
#'
#' Multiple Triangular Flow Direction:
#'
#' Seibert, J., McGlynn, B. (2007): A new triangular multiple flow direction algorithm for computing upslope areas from gridded digital elevation models. Water Ressources Research, 43, W04501.
#'
#' @author Alexander Brenning (R interface), Olaf Conrad (SAGA module), Thomas Grabs (MTFD algorithm)
#' @note This function uses module \code{Parallel Processing} (version 2.0.7+: \code{Catchment Area (Parallel)} from SAGA library \code{ta_hydrology}.
#'
#' The SAGA GIS 2.0.6+ version of the module adds more (optional) input and 
#' output grids that are currently not supported by this wrapper function.
#' Use \code{\link{rsaga.geoprocessor}} for access to these options,
#' and see \code{rsaga.get.usage("ta_hydrology","Catchment Area (Parallel)")}
#' for information on new arguments.
#' @seealso \code{\link{rsaga.topdown.processing}}, \code{\link{rsaga.wetness.index}}, \code{\link{rsaga.geoprocessor}}, \code{\link{rsaga.env}}
#' @examples
#' \dontrun{
#' # SAGA GIS 2.0.6+:
#' rsaga.get.usage("ta_hydrology","Catchment Area (Parallel)")
#' # earlier versions of SAGA GIS:
#' #rsaga.get.usage("ta_hydrology","Parallel Processing")
#' # execute model with typical settings:
#' rsaga.parallel.processing(in.dem = "dem", out.carea = "carea", out.cslope = "cslope")
#' # cslope is in radians - convert to degree:
#' fac = round(180/pi, 4)
#' formula = paste(fac, "*a", sep = "")
#' rsaga.grid.calculus("cslope", "cslopedeg", formula)
#' }
#' @keywords spatial interface
#' @export
rsaga.parallel.processing = function(in.dem, in.sinkroute, in.weight,
    out.carea, out.cheight, out.cslope, out.caspect, out.flowpath,
    step, method="mfd", linear.threshold=Inf, convergence=1.1,
    env = rsaga.env(), ...)
{
    ## Version Stop - tool no longer supported SAGA 2.1.3
    if (env$version == "2.1.3" | env$version == "2.1.4" | env$version == "2.2.0" | env$version == "2.2.1" |
        env$version == "2.2.2" | env$version == "2.2.3") {
      stop("Parallel processing not supported with SAGA GIS 2.1.3 and higher;\n",
           "See help(rsaga.topdown.processing) for similar function with SAGA 2.1.3+")  
    }
    in.dem = default.file.extension(in.dem,".sgrd")
    pp.choices = c("d8","rho8","braunschweig","dinf","mfd", "mtfd")
    method = match.arg.ext(method, choices=pp.choices,
        numeric=TRUE, ignore.case=TRUE, base=0)
    param = list( ELEVATION=in.dem )
    if (!missing(in.sinkroute)) {
        in.sinkroute = default.file.extension(in.sinkroute,".sgrd")
        param = c(param, SINKROUTE=in.sinkroute)
    }
    if (!missing(in.weight)) {
        in.weight = default.file.extension(in.weight,".sgrd")
        param = c(param, SINKROUTE=in.weight)
    }
    if (!missing(out.carea))
        param = c(param, CAREA=out.carea)
    if (!missing(out.cheight))
        param = c(param, CHEIGHT=out.cheight)
    if (!missing(out.cslope))  
        param = c(param, CSLOPE=out.cslope)
    if (!missing(step))
        param = c(param, STEP=step)
    if (!missing(out.caspect))
        param = c(param, CASPECT=out.caspect)
    if (!missing(out.flowpath))
        param = c(param, FLWPATH=out.flowpath)
    param = c(param, Method=method)
    if (is.finite(linear.threshold)) {
        param = c(param, DOLINEAR=TRUE, LINEARTHRS=linear.threshold)
    } else param = c(param, DOLINEAR=FALSE)
    
    param = c(param, CONVERGENCE=convergence)
    
    module = "Catchment Area (Parallel)"
    if (env$version == "2.0.4" | env$version == "2.0.5" | env$version == "2.0.6")
        module = "Parallel Processing"
    
    rsaga.geoprocessor(lib = "ta_hydrology", module = module, param, env = env, ...)
}

#' Top-Down Processing
#' 
#' Calculate the size of the local catchment area (contributing area), accumulated material, and flow path length, using top-down processing algorithms from the highest to the lowest cell. \cr Top-Down Processing is new with SAGA GIS 2.1.3. See \code{\link{rsaga.parallel.processing}} with older versions.
#' @name rsaga.topdown.processing
#' @param in.dem input: digital elevation model (DEM) as SAGA grid file (default file extension: \code{.sgrd})
#' @param in.sinkroute optional input: SAGA grid with sink routes
#' @param in.weight optional input: SAGA grid with weights
#' @param in.mean optional input: SAGA grid for mean over catchment calculation
#' @param in.material optional input: SAGA grid with material
#' @param in.target optional input: SAGA grid of accumulation target
#' @param in.lin.val optional input: SAGA grid providing values to be compared with linear flow threshold instead of catchment area
#' @param in.lin.dir optional input: SAGA grid to be used for linear flow routing, if the value is a valid direction (0-7 = N, NE, E, SE, S, SW, W, NW)
#' @param out.carea output: catchment area grid
#' @param out.mean optional output: mean over catchment grid
#' @param out.tot.mat optional output: total accumulated material grid
#' @param out.acc.left optional output: accumulated material from left side grid
#' @param out.acc.right optional output: accumulated material from right side grid
#' @param out.flowpath optional output: flow path length grid
#' @param step integer >=1: step parameter
#' @param method character or numeric: choice of processing algorithm (default \code{"mfd"}, or 4):
#' \itemize{
#' \item [0] Deterministic 8 (\code{"d8"} or 0)
#' \item [1] Rho 8 (\code{"rho8"}, or 1)
#' \item [2] Braunschweiger Reliefmodell (\code{"braunschweig"} or 2)
#' \item [3] Deterministic Infinity (\code{"dinf"} or 3)
#' \item [4] Multiple Flow Direction (\code{"mfd"} or 4)
#' \item [5] Multiple Triangular Flow Direction (\code{"mtfd"}, or 5)
#' \item [6] Multiple Maximum Gradient Based Flow Direction (\code{"mdg"}, or 6)}
#' @param linear.threshold numeric (number of grid cells): threshold above which linear flow (i.e. the Deterministic 8 algorithm) will be used; linear flow is disabled for \code{linear.threshold=Inf} (the default)
#' @param convergence numeric >=0: a parameter for tuning convergent/ divergent flow; default value of \code{1.1} gives realistic results and should not be changed
#' @param env list, setting up a SAGA geoprocessing environment as created by \code{\link{rsaga.env}}
#' @param ... further arguments to \code{\link{rsaga.geoprocessor}}
#' @details Refer to the references for details on the available algorithms.
#' @return The type of object returned depends on the \code{intern} argument passed to the \code{\link{rsaga.geoprocessor}}. For \code{intern=FALSE} it is a numerical error code (0: success), or otherwise (the default) a character vector with the module's console output.
#' @references
#' Deterministic 8:
#'
#' O'Callaghan, J.F., Mark, D.M. (1984): The extraction of drainage networks from digital elevation data. Computer Vision, Graphics and Image Processing, 28: 323-344.
#'
#' Rho 8:
#'
#' Fairfield, J., Leymarie, P. (1991): Drainage networks from grid digital elevation models. Water Resources Research, 27: 709-717.
#'
#' Braunschweiger Reliefmodell:
#'
#' Bauer, J., Rohdenburg, H., Bork, H.-R. (1985): Ein Digitales Reliefmodell als Vorraussetzung fuer ein deterministisches Modell der Wasser- und Stoff-Fluesse. Landschaftsgenese und Landschaftsoekologie, H. 10, Parameteraufbereitung fuer deterministische Gebiets-Wassermodelle, Grundlagenarbeiten zu Analyse von Agrar-Oekosystemen, eds.: Bork, H.-R., Rohdenburg, H., p. 1-15.
#'
#' Deterministic Infinity:
#'
#' Tarboton, D.G. (1997): A new method for the determination of flow directions and upslope areas in grid digital elevation models. Water Ressources Research, 33(2): 309-319.
#'
#' Multiple Flow Direction:
#'
#' Freeman, G.T. (1991): Calculating catchment area with divergent flow based on a regular grid. Computers and Geosciences, 17: 413-22.
#'
#' Quinn, P.F., Beven, K.J., Chevallier, P., Planchon, O. (1991): The prediction of hillslope flow paths for distributed hydrological modelling using digital terrain models. Hydrological Processes, 5: 59-79.
#'
#' Multiple Triangular Flow Direction:
#'
#' Seibert, J., McGlynn, B. (2007): A new triangular multiple flow direction algorithm for computing upslope areas from gridded digital elevation models. Water Ressources Research, 43, W04501.
#'
#' Multiple Flow Direction Based on Maximum Downslope Gradient:
#' 
#' Qin, C.Z., Zhu, A-X., Pei, T., Li, B.L., Scholten, T., Zhou, C.H. (2011): An approach to computing topographic wetness index based on maximum downslope gradient. Precision Agriculture, 12(1): 32-43.
#' 
#' @author Alexander Brenning and Donovan Bangs (R interface), Olaf Conrad (SAGA module), Thomas Grabs (MTFD algorithm)
#' @examples
#' \dontrun{
#' # Calculation of contributing area with default settings:
#' rsaga.topdown.processing(in.dem = "dem", out.carea = "carea")
#' # Calculation of contributing area by maximunm downslope gradient:
#' rsaga.topdown.processing(in.dem = "dem", out.carea = "carea",
#'                          method = "mdg")
#' }
#' @seealso \code{\link{rsaga.parallel.processing}}, \code{\link{rsaga.wetness.index}}, \code{\link{rsaga.geoprocessor}}, \code{\link{rsaga.env}}
#' @keywords spatial interface
#' @export
rsaga.topdown.processing = function(in.dem, in.sinkroute, in.weight, in.mean, in.material, in.target,
                                    in.lin.val, in.lin.dir,
                                    out.carea, out.mean, out.tot.mat, out.acc.left, out.acc.right,
                                    out.flowpath, step, method = "mfd", linear.threshold = Inf, convergence = 1.1,
                                    env = rsaga.env(), ...) {
    ## Version Stop - SAGA GIS Version < 2.1.3
    if (env$version != "2.1.3" & env$version != "2.1.4" & env$version != "2.2.0" & env$version != "2.2.1" &
        env$version != "2.2.2" & env$version != "2.2.3") {
        stop("rsaga.topdown.processing requires SAGA GIS 2.1.3 or higher;\n",
             "see help(rsaga.parallel.processing) for similar function in earlier versions")
    }
    
    in.dem = default.file.extension(in.dem,".sgrd")
    pp.choices = c("d8","rho8","braunschweig","dinf","mfd", "mtfd", "mdg")
    method = match.arg.ext(method, choices=pp.choices,
                           numeric=TRUE, ignore.case=TRUE, base=0)
    param = list( ELEVATION=in.dem )
    if (!missing(in.sinkroute)) {
        in.sinkroute = default.file.extension(in.sinkroute,".sgrd")
        param = c(param, SINKROUTE=in.sinkroute)
    }
    if (!missing(in.weight)) {
        in.weight = default.file.extension(in.weight,".sgrd")
        param = c(param, SINKROUTE=in.weight)
    }
    if (!missing(in.mean)) {
        in.mean = default.file.extension(in.mean, ".sgrd")
        param = c(param,VAL_INPUT=in.mean)
    }
    if (!missing(in.material)) {
        in.material = default.file.extension(in.material, ".sgrd")
        param = c(param, MATERIAL=in.material)
    }
    if (!missing(in.target)) {
        in.target = default.file.extension(in.target, ".sgrd")
        param = c(param, TARGET=in.target)
    }
    if (!missing(in.lin.val)) {
        in.lin.val = default.file.extension(in.lin.val, ".sgrd")
        param = c(param, LINEAR_VAL=in.lin.val)
    }
    if (!missing(in.lin.dir)){
        in.lin.dir = default.file.extension(in.lin.dir, ".sgrd")
        param = c(param, LINEAR_DIR=in.lin.dir)
    }
    if (!missing(out.carea))
        param = c(param, CAREA=out.carea)
    if (!missing(out.mean))
        param = c(param, VAL_MEAN=out.mean)
    if (!missing(out.tot.mat))
        param = c(param, ACCU_TOT=out.tot.mat)
    if (!missing(out.acc.left))
        param = c(param, ACCU_LEFT=out.acc.left)
    if (!missing(out.acc.right))
        param = c(param, ACCU_RIGHT=out.acc.right)
    if (!missing(out.flowpath))
        param = c(param, FLOWLEN=out.flowpath)
    param = c(param, METHOD=method)
    if (is.finite(linear.threshold)) {
        param = c(param, LINEAR_DO=TRUE, LINEAR_MIN=linear.threshold)
    } else param = c(param, LINEAR_DO=FALSE)
    
    param = c(param, CONVERGENCE=convergence)
    
    module = "Catchment Area (Top-Down)"

    if (env$version == "2.2.0" | env$version == "2.2.1" | env$version == "2.2.2" |
        env$version == "2.2.3") {
        module = "Flow Accumulation (Top-Down)"
    }

    rsaga.geoprocessor(lib = "ta_hydrology", module = module, param, env = env, ...) 
}

#' SAGA Modules SAGA Wetness Index
#' 
#' Calculate the SAGA Wetness Index (SWI), a modified topographic wetness index (TWI)
#' @name rsaga.wetness.index
#' @param in.dem input: digital elevation model (DEM) as SAGA grid file (default file extension: \code{.sgrd})
#' @param out.wetness.index output file (optional): wetness index grid file name. Existing files of the same name will be overwritten!
#' @param out.carea output file (optional): catchment area grid file name
#' @param out.cslope output file (optional): catchment slope grid file name
#' @param out.mod.carea output file (optional): file name of modified catchment area grid
#' @param suction SAGA GIS 2.1.0+: positive numeric value (optional): the lower this value is the stronger is the suction effect; defaults to a value of 10 (more detailed information is currently not available  in the SAGA GIS documentation
#' @param area.type character or numeric (optional): type of area: \code{"absolute"} (or numeric code 0): absolute catchment area; \code{"square root"} (code 1; the default): square root of catchment area; \code{"specific"} (code 2): specific catchment area
#' @param slope.type character or numeric (optional): type of slope: \code{"local"} (or numeric code 0): local slope; \code{"catchment"} (or code 1; the default): catchment slope.
#' @param slope.min numeric (optional): minimum slope; default: 0
#' @param slope.offset numeric (optional): offset slope; default: 0.1
#' @param slope.weight numeric (optional): weighting factor for slope in index calculation; default: 1
#' @param t.param SAGA GIS up to version 2.0.8: positive numeric value (optional): undocumented
#' @param env A SAGA geoprocessing environment, see \code{\link{rsaga.env}}.)
#' @param ... optional arguments to be passed to \code{\link{rsaga.geoprocessor}} 
#' @details The SAGA Wetness Index is similar to the  Topographic Wetness Index (TWI), but it is based on a modified  catchment area calculation (\code{out.mod.carea}), which does not treat the flow as a thin film as done in the calculation of catchment areas in conventional algorithms. As a result, the SWI tends to assign a more realistic, higher potential soil wetness than the TWI to grid cells situated in valley floors with a small vertical distance to a channel.
#'
#' This module and its arguments changed substantially from SAGA GIS 2.0.8 to version 2.1.0. It appears to me that the new algorithm is similar (but not identical) to the old one when using \code{area.type="absolute"} and \code{slope.type="local"} but I haven't tried out all possible options. This help file will be updated as soon as additional documentation becomes available.
#' @return The type of object returned depends on the \code{intern} argument passed to the \code{\link{rsaga.geoprocessor}}. For \code{intern=FALSE} it is a numerical error code (0: success), or otherwise (the default) a character vector with the module's console output.
#' @references Boehner, J., Koethe, R. Conrad, O., Gross, J.,  Ringeler, A., Selige, T. (2002): Soil Regionalisation by Means of Terrain Analysis and Process Parameterisation. In: Micheli, E., Nachtergaele, F., Montanarella, L. (ed.): Soil Classification 2001. European Soil Bureau, Research Report No. 7, EUR 20398 EN, Luxembourg. pp.213-222.
#' 
#' Boehner, J. and Selige, T. (2006): Spatial prediction of soil attributes using terrain analysis and climate regionalisation. In: Boehner, J., McCloy, K.R., Strobl, J. [Ed.]: SAGA - Analysis and Modelling Applications, Goettinger Geographische Abhandlungen, Goettingen: 13-28.
#' @author Alexander Brenning (R interface), Juergen Boehner and Olaf Conrad (SAGA module)
#' @seealso \code{\link{rsaga.parallel.processing}}, \code{\link{rsaga.geoprocessor}}, \code{\link{rsaga.env}}
#' @examples
#' \dontrun{
#' # using SAGA grids:
#' rsaga.wetness.index("dem.sgrd","swi.sgrd")
#' }
#' @keywords spatial interface
#' @export
rsaga.wetness.index = function( in.dem, 
    out.wetness.index, out.carea, out.cslope, 
    out.mod.carea, 
    # since SAGA GIS 2.1.0:
    suction, area.type, slope.type, slope.min, slope.offset, slope.weight,
    # up to SAGA GIS 2.0.8:
    t.param,
    env = rsaga.env(), ...)
{
    in.dem = default.file.extension(in.dem,".sgrd")
    if (missing(out.carea)) {
        out.carea = tempfile()
        on.exit(unlink(paste(out.carea,".*",sep="")), add = TRUE)
    }
    if (missing(out.cslope)) {
        out.cslope = tempfile()
        on.exit(unlink(paste(out.cslope,".*",sep="")), add=TRUE)
    }
    if (missing(out.mod.carea)) {
        out.mod.carea = tempfile()
        on.exit(unlink(paste(out.mod.carea,".*",sep="")), add=TRUE)
    }
    if (env$version == "2.1.0" | env$version == "2.1.1" | env$version == "2.1.2" |
        env$version == "2.1.3" | env$version == "2.1.4" | env$version == "2.2.0" |
        env$version == "2.2.1" | env$version == "2.2.2" | env$version == "2.2.3")  {
        param = list(DEM=in.dem, AREA=out.carea, SLOPE=out.cslope, 
                     AREA_MOD=out.mod.carea, TWI=out.wetness.index)
        if (!missing(suction)) {
            suction = as.numeric(suction)
            if (suction <= 0) stop("'suction' argument must be >0")
            param = c(param, SUCTION=suction)
        }
        if (!missing(area.type)) {
            area.type = match.arg.ext(area.type,choices=c("absolute","square root","specific"),base=0,ignore.case=TRUE,numeric=TRUE)
            param = c(param, AREA_TYPE=area.type)
        }
        if (!missing(slope.type)) {
            slope.type = match.arg.ext(slope.type,choices=c("local","catchment"),base=0,ignore.case=TRUE,numeric=TRUE)
            param = c(param, SLOPE_TYPE=slope.type)
        }
        if (!missing(slope.min)) {
            slope.min = as.numeric(slope.min)
            if (slope.min < 0) stop("'slope.min' argument must be >=0")
            param = c(param, SLOPE.MIN=slope.min)
        }
        if (!missing(slope.offset)) {
            slope.offset = as.numeric(slope.offset)
            if (slope.offset < 0) stop("'slope.offset' argument must be >=0")
            param = c(param, SLOPE.OFF=slope.offset)
        }
        if (!missing(slope.weight)) {
            slope.weight = as.numeric(slope.weight)
            if (slope.weight < 0) stop("'slope.weight' argument must be >=0")
            param = c(param, SLOPE.WEIGHT=slope.weight)
        }
        if (!missing(t.param))
            warning("argument 't.param' (in saga_cmd: T) supported only up to SAGA GIS 2.0.8")
    } else {
        param = list(DEM=in.dem, C=out.carea, GN=out.cslope, 
                     CS=out.mod.carea, SB=out.wetness.index)
        if (!missing(t.param))
            param = c(param, T=as.numeric(t.param))
        if (!missing(suction) | !missing(area.type) | !missing(slope.type) | !missing(slope.min) | !missing(slope.offset) | !missing(slope.weight))
            warning("arguments 'suction', 'area.type', 'slope.min', 'slope.type', 'slope.offset'\n",
                    "and 'slope.weight' not supported prior to SAGA GIS 2.1.0")
    }
    rsaga.geoprocessor(lib = "ta_hydrology",
        module = "SAGA Wetness Index",
        param, ..., env = env)
}






########    Module grid_calculus    ########



#' SAGA Module Grid Calculus
#'
#' Perform Arithmetic Operations on Grids
#' @name rsaga.grid.calculus
#' @param in.grids input character vector: SAGA grid files (default file extension: \code{.sgrd})
#' @param out.grid output: grid file resulting from the cell-by-cell application of 'formula' to the grids. Existing files will be overwritten!
#' @param formula character string of formula specifying the arithmetic operation to be performed on the \code{in.grids} (see Details); if this is a formula, only the right hand side will be used.
#' @param coef numeric: coefficient vector to be used for the linear combination of the \code{in.grids}. If \code{coef} as one more element than \code{in.grids}, the first one will be interpreted as an intercept.
#' @param cf.digits integer: number of digits used when converting the \code{coef}ficients to character strings (trailing zeros will be removed)
#' @param remove.zeros logical: if \code{TRUE}, terms (grids) with coefficient (numerically) equal to zero (after rounding to \code{cf.digits} digits) will be removed from the formula
#' @param remove.ones logical: if \code{TRUE} (th edefault), factors equal to 1 (after rounding to \code{cf.digits} digits) will be removed from the formula
#' @param env RSAGA geoprocessing environment, generated by a call to \code{\link{rsaga.env}}
#' @param ... optional arguments to be passed to \code{\link{rsaga.geoprocessor}}
#' @details The \code{in.grids} are represented in the \code{formula} by the letters \code{a} (for \code{in.grids[1]}), \code{b} etc. Thus, if \code{in.grids[1]} is Landsat TM channel 3 and \code{in.grids[2]} is channel 4, the NDVI formula (TM3-TM4)/(TM3+TM4) can be represented  by the character string \code{"(a-b)/(a+b)"} (any spaces are removed) or the formula \code{~(a-b)/(a+b)} in the \code{formula} argument.
#' 
#' In addition to +, -, *, and /, the following operators and functions are available for the \code{formula} definition:
#' \itemize{
#'     \item \eqn{\hat{\ }}{^} power
#'     \item \code{sin(a)} sine
#'     \item \code{cos(a)} cosine
#'     \item \code{tan(a)} tangent
#'     \item \code{asin(a)} arc sine
#'     \item \code{acos(a)} arc cosine
#'     \item \code{atan(a)} arc tangent
#'     \item \code{atan2(a,b)} arc tangent of b/a
#'     \item \code{abs(a)} absolute value
#'     \item \code{int(a)} convert to integer
#'     \item \code{sqr(a)} square
#'     \item \code{sqrt(a)} square root
#'     \item \code{ln(a)} natural logarithm
#'     \item \code{log(a)} base 10 logarithm
#'     \item \code{mod(a,b)} modulo
#'     \item \code{gt(a, b)} returns 1 if a greater b
#'     \item \code{lt(a, b)} returns 1 if a lower b
#'     \item \code{eq(a, b)} returns 1 if a equal b
#'     \item \code{ifelse(switch, x, y)} returns x if switch equals 1 else y
#' }
#' 
#' Using \code{remove.zeros=FALSE} might have the side effect that no data areas in the grid with coefficient 0 are passed on to the results grid. (To be confirmed.)
#' @return The type of object returned depends on the \code{intern} argument passed to the \code{\link{rsaga.geoprocessor}}. For \code{intern=FALSE} it is a numerical error code (0: success), or otherwise (the default) a character vector with the module's console output.
#' @author Alexander Brenning (R interface), Olaf Conrad (SAGA module)
#' @seealso \code{\link{local.function}}, \code{\link{focal.function}}, and \code{\link{multi.focal.function}} for a more flexible framework for combining grids or applying local and focal functions; \code{\link{rsaga.geoprocessor}}, \code{\link{rsaga.env}}
#' @examples
#' \dontrun{
#' # using SAGA grids:
#' # calculate the NDVI from Landsat TM bands 3 and 4:
#' rsaga.grid.calculus(c("tm3.sgrd","tm4.sgrd"), "ndvi.sgrd", ~(a-b)/(a+b))
#' # apply a linear regression equation to grids:
#' coefs = c(20,-0.6)
#' # maybe from a linear regression of mean annual air temperature (MAAT)
#' # against elevation - something like:
#' # coefs = coef( lm( maat ~ elevation ) )
#' rsaga.linear.combination("elevation.sgrd", "maat.sgrd", coefs)
#' # equivalent:
#' rsaga.grid.calculus("elevation.sgrd", "maat.sgrd", "20 - 0.6*a")
#' }
#' @keywords spatial interface
#' @export
rsaga.grid.calculus = function(in.grids, out.grid, formula,
    env = rsaga.env(), ...)
{
    in.grids = default.file.extension(in.grids, ".sgrd")
    in.grids = paste(in.grids, collapse = ";")
    if (any(class(formula) == "formula"))
        formula = rev( as.character(formula) )[1]
    formula = gsub(" ", "", formula)
    if (env$version == "2.0.4") {
        param = list( INPUT = in.grids, RESULT = out.grid,
                    FORMUL = formula )
    } else {
        param = list( GRIDS = in.grids, RESULT = out.grid,
                    FORMULA = formula )
    }
    rsaga.geoprocessor(lib = "grid_calculus", 
        module = "Grid Calculator", # was = 1
        param = param, env = env, ...)
}



#' @rdname rsaga.grid.calculus
#' @name rsaga.linear.combination
#' @export
rsaga.linear.combination = function(in.grids, out.grid, coef, 
    cf.digits = 16, remove.zeros = FALSE, remove.ones = TRUE, 
    env = rsaga.env(), ...)
{
    fmt = paste("%.", cf.digits, "f", sep = "")
    coef = sprintf(fmt, coef)
    zero = sprintf(fmt, 0)
    omit = rep(FALSE, length(coef))

    if (length(coef) == length(in.grids)) { # no intercept provided
        coef = c(NA, coef)
        omit = c(TRUE, omit)
    }
    nvars = length(coef)
    if (nvars != length(in.grids) + 1)
        stop("'coef' must have length 'length(in.grids)' or 'length(in.grids)+1'")

    # Simplify the formula by removing terms that are zero
    # (after rounding to the specified number of digits):
    if (remove.zeros)
        omit = omit | (coef == zero)
    # Zero intercept is always removed:
    omit[1] = omit[1] | (coef[1] == zero)

    # Remove zeros at the end of the coefficients:
    for (i in 1:nvars) {
        if (omit[i]) next
        # Are there any digits at all?
        if (length(grep(".", coef[i], fixed = TRUE)) == 0) next
        nc = nchar(coef[i])
        # Remove all trailing zeros:
        while (substr(coef[i], nc, nc) == "0") {
            coef[i] = substr(coef[i], 1, nc - 1)
            nc = nchar(coef[i])
        }
        # Remove trailing decimal point:
        if (substr(coef[i], nc, nc) == ".")
            coef[i] = substr(coef[i], 1, nc - 1)
    }

    # Set up the formula:
    ltrs = letters[ 1 : sum(!omit[-1]) ]
    if (!omit[1]) ltrs = c("intercept", ltrs)
    formula = paste(coef[ !omit ], ltrs, 
                    collapse = "+", sep = "*")
    formula = gsub("*intercept", "", formula, fixed = TRUE)
    formula = gsub("+-", "-", formula, fixed = TRUE)
    if (remove.ones) {
        formula = gsub("-1*", "-", formula, fixed = TRUE)
        formula = gsub("+1*", "+", formula, fixed = TRUE)
    }
    
    rsaga.grid.calculus(in.grids = in.grids[!omit[-1]], out.grid = out.grid,
        formula = formula, env = env, ...)
}





########     Module shapes_grid     ########



#' Contour Lines from a Grid
#'
#' Creates a contour lines shapefile from a grid file in SAGA grid format.
#' @name rsaga.contour
#' @param in.grid input: digital elevation model (DEM) as SAGA grid file (default file extension: \code{.sgrd})
#' @param out.shapefile output: contour line shapefile. Existing files will be overwritten!
#' @param zstep,zmin,zmax lower limit, upper limit, and equidistance of contour lines
#' @param vertex optional parameter: vertex type for resulting contours. Default \code{"xy"} (or 0). Only available with SAGA GIS 2.1.3+. \itemize{
#' \item [0] \code{"xy"}
#' \item [1] \code{"xyz"}}
#' @param env A SAGA geoprocessing environment, see \code{\link{rsaga.env}}
#' @param ... arguments to be passed to \code{\link{rsaga.geoprocessor}}
#' @return The type of object returned depends on the \code{intern} argument passed to the \code{\link{rsaga.geoprocessor}}. For \code{intern=FALSE} it is a numerical error code (0: success), or otherwise (the default) a character vector with the module's console output.
#' @author Alexander Brenning (R interface), Olaf Conrad (SAGA module)
#' @seealso \code{\link{rsaga.geoprocessor}}
#' @keywords spatial interface
#' @export
rsaga.contour = function(in.grid,out.shapefile,zstep,zmin,zmax,vertex="xy",env=rsaga.env(),...) {
    in.grid = default.file.extension(in.grid,".sgrd")
    # 'INPUT' changed to 'GRID' with SAGA 2.1.3
    if(env$version != "2.1.3" & env$version != "2.1.4" & env$version != "2.2.0" & env$version != "2.2.1" &
       env$version != "2.2.2" & env$version != "2.2.3"){
        param = list(INPUT=in.grid,CONTOUR=out.shapefile)
    } else {
        param = list(GRID=in.grid,CONTOUR=out.shapefile)
    }
    if (!missing(zmin))  param = c(param, ZMIN=as.numeric(zmin))
    if (!missing(zmax))  param = c(param, ZMAX=as.numeric(zmax))
    if (!missing(zstep)) {
        stopifnot(as.numeric(zstep)>0)
        param = c(param, ZSTEP=as.numeric(zstep))
    }
    v.choices = c("xy", "xyz")
    vertex = match.arg.ext(vertex, choices=v.choices,
                           numeric=TRUE, ignore.case=TRUE, base=0)
    if (!missing(vertex)) {
        if (env$version == "2.1.3" | env$version == "2.1.4") {
            param = c(param, VERTEX=vertex)
        }
    }
    rsaga.geoprocessor(lib = "shapes_grid", 
        module = "Contour Lines from Grid",
        param, env = env,...)
}



#' Add Grid Values to Point Shapefile
#'
#' Pick values from SAGA grids and attach them as a new variables to a point shapefile.
#' @name rsaga.add.grid.values.to.points
#' @param in.grids Input: character vector with names of (one or more) SAGA GIS grid files to be converted into a point shapefile.
#' @param in.shapefile Input point shapefile (default extension: \code{.shp}).
#' @param out.shapefile Output point shapefile (default extension: \code{.shp}).
#' @param method interpolation method to be used; choices: nearest neighbour interpolation (default), bilinear interpolation, inverse distance weighting, bicubic spline interpolation, B-splines.
#' @param ... Optional arguments to be passed to \code{\link{rsaga.geoprocessor}}, including the \code{env} RSAGA geoprocessing environment.
#' @details Retrieves information from the selected grids at the positions of the points of the selected points layer and adds it to the resulting layer.
#' @author Alexander Brenning (R interface), Olaf Conrad (SAGA modules)
#' @note This function uses module \code{Add Grid Values to Points} in SAGA GIS library \code{shapes_grid}.
#' @seealso \code{\link{pick.from.points}}, \code{\link{pick.from.ascii.grid}}, \code{\link{pick.from.saga.grid}}, \code{\link{rsaga.grid.to.points}}
#' @keywords spatial interface
#' @export
rsaga.add.grid.values.to.points = function(in.shapefile,
    in.grids, out.shapefile, 
    method = c("nearest.neighbour", "bilinear",
      "idw", "bicubic.spline", "b.spline"), ...)
{
    in.grids = default.file.extension(in.grids,".sgrd")
    in.grids = paste(in.grids, collapse = ";")
    # check if this is SAGA version dependent:
    in.shapefile = default.file.extension(in.shapefile,".shp")
    out.shapefile = default.file.extension(out.shapefile,".shp")
    method = match.arg.ext(method, base = 0, ignore.case = TRUE, numeric = TRUE)
    param = list(SHAPES = in.shapefile, GRIDS = in.grids,
                RESULT = out.shapefile, INTERPOL = method)
    rsaga.geoprocessor(lib = "shapes_grid", 
        module = "Add Grid Values to Points", # was: = 0
        param, ...)
}


#' Convert SAGA grid file to point shapefile
#' 
#' Convert SAGA grid file to point (or polygon) shapefile - either completely or only a random sample of grid cells.
#' @name rsaga.grid.to.points
#' @param in.grids Input: names of (possibly several) SAGA GIS grid files to be converted into a point shapefile.
#' @param in.grid Input: SAGA grid file from which to sample.
#' @param out.shapefile Output: point shapefile (default extension: \code{.shp}). Existing files will be overwritten!
#' @param in.clip.polygons optional polygon shapefile to be used for clipping/masking an area
#' @param exclude.nodata logical (default: \code{TRUE}): skip 'nodata' grid cells?
#' @param type character string: \code{"nodes"}: create point shapefile of grid center points; \code{"cells"} (only supported by SAGA GIS 2.0.6+): create polygon shapefile with grid cell boundaries
#' @param freq integer >=1: sampling frequency: on average 1 out of 'freq' grid cells are selected
#' @param env RSAGA geoprocessing environment created by \code{\link{rsaga.env}}; required by \code{rsaga.grid.to.points} to determine version-dependent SAGA module name and arguments
#' @param ... Optional arguments to be passed to \code{\link{rsaga.geoprocessor}}
#' @author Alexander Brenning (R interface), Olaf Conrad (SAGA modules)
#' @note These functions use modules \code{Grid Values to Shapes} (pre-2.0.6 name: \code{Grid Values to Points}) and \code{Grid Values to Points (randomly)} in SAGA library \code{shapes_grid}.
#'
#' The SAGA 2.0.6+ module \code{Grid Values to Shapes} is more flexible than the earlier versions as it allows to create grid cell polygons instead of center points (see argument \code{type}).
#' @seealso \code{\link{rsaga.add.grid.values.to.points}}
#' @examples
#' \dontrun{
#' # one point per grid cell, exclude nodata areas:
#' rsaga.grid.to.points("dem", "dempoints")
#' # take only every 20th point, but to not exclude nodata areas:
#' rsaga.grid.to.points.randomly("dem", "dempoints20", freq = 20)
#' }
#' @keywords spatial interface
#' @export
rsaga.grid.to.points = function(in.grids, out.shapefile, 
    in.clip.polygons, exclude.nodata = TRUE,
    type = "nodes", env = rsaga.env(), ...)
{
    in.grids = default.file.extension(in.grids,".sgrd")
    in.grids = paste(in.grids, collapse = ";")
    type = match.arg.ext(type, numeric=TRUE, ignore.case=TRUE, base=0,
        choices=c("nodes","cells"))
    if (type == 1 & (env$version == "2.0.4" | env$version == "2.0.5")) {
        type = 0
        warning("type == 'cells' not supported by SAGA 2.0.4 and 2.0.5; using type = 'nodes'")
    }
    param = list(GRIDS = in.grids)
    if (env$version == "2.0.4" | env$version == "2.0.5") {
        param = c(param, POINTS = out.shapefile)
    } else param = c(param, SHAPES = out.shapefile)
    param = c(param, NODATA = exclude.nodata)
    if (!missing(in.clip.polygons))
        param = c(param, POLYGONS = in.clip.polygons)
    if (!(env$version == "2.0.4" | env$version == "2.0.5"))
        param = c(param, TYPE = type)
    module = "Grid Values to Shapes"
    if (!rsaga.module.exists("shapes_grid",module,env=env))
    #if (env$version == "2.0.4" | env$version == "2.0.5")
        module = "Grid Values to Points"
    rsaga.geoprocessor(lib = "shapes_grid", 
        module = module, # was: = 3
        param, env = env, ...)
}


#' @rdname rsaga.grid.to.points
#' @name rsaga.grid.to.points.randomly
#' @export
rsaga.grid.to.points.randomly = function(in.grid,
    out.shapefile, freq, ...)
{
    in.grid = default.file.extension(in.grid, ".sgrd")
    out.shapefile = default.file.extension(out.shapefile, ".shp")
    if (freq < 1) stop("'freq' must be an integer >=1")
    param = list(GRID = in.grid, FREQ = freq, POINTS = out.shapefile)
    rsaga.geoprocessor(lib = "shapes_grid", 
        module = "Grid Values to Points (randomly)", # was: = 4
        param, ...)
}



#' Spatial Interpolation Methods
#' 
#' Spatial interpolation of point data using inverse distance to a power (inverse distance weighting, IDW), nearest neighbors, or modified quadratic shephard.
#' @name rsaga.inverse.distance
#' @param in.shapefile Input: point shapefile (default extension: \code{.shp}).
#' @param out.grid Output: filename for interpolated grid (SAGA grid file). Existing files will be overwritten!
#' @param field numeric or character: number or name of attribute in the shapefile's attribute table to be interpolated; the first attribute is represented by a zero.
#' @param power numeric (>0): exponent used in inverse distance  weighting (usually 1 or 2)
#' @param maxdist numeric: maximum distance of points to be used for inverse distance interpolation (search radius); no search radius is applied when this argument is missing or equals \code{Inf}
#' @param nmax Maximum number of nearest points to be used for interpolation; \code{nmax=Inf} is a valid value (no upper limit)
#' @param quadratic.neighbors integer >=5; default 13.
#' @param weighting.neighbors integer >=3; default 19.
#' @param target required argument of type list: parameters identifying the target area, e.g. the x/y extent and cellsize, or name of a reference grid; see \code{\link{rsaga.target}}.
#' @param env RSAGA geoprocessing environment created by \code{\link{rsaga.env}}, required because module(s) depend(s) on SAGA version
#' @param ... Optional arguments to be passed to \code{\link{rsaga.geoprocessor}}, including the \code{env} RSAGA geoprocessing environment.
#' @details These functions use modules from the \code{grid_gridding} SAGA GIS library. They do not support SAGA GIS 2.0.4, which differs in some argument names and parameterizations. Target grid parameterization by grid file name currently doesn't work with SAGA GIS 2.1.0  Release Candidate 1 (see also \code{\link{rsaga.target}}); stay tuned for future updates and fixes.
#' @references QSHEP2D: Fortran routines implementing the Quadratic Shepard method for bivariate interpolation of scattered data  (see R. J. Renka, ACM TOMS 14 (1988) pp.149-150). Classes: E2b. Interpolation of scattered, non-gridded  multivariate data.
#' @author Alexander Brenning (R interface), Andre Ringeler and Olaf Conrad (SAGA modules)
#' @note The 'Inverse Distance Weighted' module of SAGA GIS not only support inverse-distance weighted interpolation, but also exponential and other weighting schemes (command line argument WEIGHTING); these are however not accessible through this function, but only through the \code{rsaga.geoprocessor}, if needed. See \code{rsaga.get.usage("grid_gridding","Inverse Distance Weighted")} for details.
#'
#' See the example section in the help file for \code{\link[shapefiles]{write.shapefile}} in package \code{shapefiles} to learn how to apply these interpolation functions to a shapefile exported from a data.frame.
#' 
#' Modified Quadratic Shephard method: based on module 660 in TOMS (see references).
#' @seealso \code{\link{rsaga.target}}; \code{\link[gstat]{idw}} in package \code{gstat}.
#' @keywords spatial interface
#' @export
rsaga.inverse.distance = function(in.shapefile, out.grid, field, 
        power = 1, maxdist, nmax = 100,
        target, env = rsaga.env(), ...)
{
    if (env$version == "2.0.4")
        stop("rsaga.inverse.distance doesn't support SAGA GIS 2.0.4 any longer\n",
             "  because some of the arguments have changed")

    stopifnot(!missing(target))

    if (power <= 0) stop("'power' must be >0")
    if (field < 0) stop("'field' must be an integer >=0")

    in.shapefile = default.file.extension(in.shapefile, ".shp")
    out.grid = default.file.extension(out.grid, ".sgrd")
    
    if (target$TARGET == 1) {
        if (target$GRID_GRID != out.grid) {
            rsaga.copy.sgrd(target$GRID_GRID, out.grid, env = env)
            target$GRID_GRID = out.grid
        }
    }

    module = "Inverse Distance Weighted"

    param = list(
        USER_GRID = out.grid,
        SHAPES = in.shapefile,
        FIELD = field,
        WEIGHTING = 0, # IDW
        MODE = 0, # search mode: all directions
        POWER = power)

    is.global = (missing(maxdist))
    if (!missing(maxdist)) {
        if (maxdist <= 0) stop("'maxdist' must be >0")
        if (maxdist == Inf) is.global = TRUE
    }
    if (is.global) {
        param = c(param, list(RANGE = 1))
    } else
        param = c(param, list(RANGE = 0, RADIUS = maxdist))

    #use.all = (missing(nmax))
    #if (!missing(nmax)) {
    if (nmax <= 0) stop("'nmax' must be an integer >0, or Inf")
    use.all = (nmax == Inf)
    #}
    if (use.all) {
        param = c(param, list(POINTS = 1))
    } else
        param = c(param, list(POINTS = 0, NPOINTS = nmax))

    param = c(param, target)        
    
    # Translate some argument names for SAGA GIS 2.1.0+:
    if (substr(env$version,1,4) != "2.0.") {
        nm = names(param)
        nm[ nm == "RANGE" ] = "SEARCH_RANGE"
        nm[ nm == "RADIUS" ] = "SEARCH_RADIUS"
        nm[ nm == "POINTS" ] = "SEARCH_POINTS_ALL"
        nm[ nm == "NPOINTS" ] = "SEARCH_POINTS_MAX"
        nm[ nm == "MODE" ] = "SEARCH_DIRECTION"
        nm[ nm == "POWER" ] = "WEIGHT_POWER"
        # TARGET parameters changed SAGA 2.1.3:
        if (env$version == "2.1.3" | env$version == "2.1.4" | env$version == "2.2.0" |
            env$version == "2.2.1" | env$version == "2.2.2" | env$version == "2.2.3") {
            nm[ nm == "USER_GRID" ] = "TARGET_OUT_GRID"
            nm[ nm == "TARGET" ] = "TARGET_DEFINITION"
            nm[ nm == "GRID_GRID" ] = "TARGET_TEMPLATE"
            nm[ nm == "USER_SIZE" ] = "TARGET_USER_SIZE"
            nm[ nm == "USER_FIT" ] = "TARGET_USER_FITS"
            nm[ nm == "USER_XMIN" ] = "TARGET_USER_XMIN"
            nm[ nm == "USER_XMAX" ] = "TARGET_USER_XMAX"
            nm[ nm == "USER_YMIN" ] = "TARGET_USER_YMIN"
            nm[ nm == "USER_YMAX" ] = "TARGET_USER_YMAX"
        }
        names(param) = nm
        
        # Translate some argument names for SAGA 2.2.0
        if (substr(env$version,1,4) == "2.2."){
            nm = names(param)
            nm[ nm == "WEIGHTING" ] = "DW_WEIGHTING"
            nm[ nm == "WEIGHT_POWER" ] = "DW_IDW_POWER"
            nm[ nm == "WEIGHT_BANDWIDTH" ] = "DW_BANDWIDTH"
        }
        
        names(param) = nm
    }

    rsaga.geoprocessor(lib = "grid_gridding", 
        module = module,
        param = param, env = env, ...)
}


#' @rdname rsaga.inverse.distance
#' @name rsaga.nearest.neighbour
#' @export
rsaga.nearest.neighbour = function(in.shapefile, out.grid, field,
    target, env = rsaga.env(), ...)
{
    if (env$version == "2.0.4")
        stop("rsaga.nearest.neighbour doesn't support SAGA GIS 2.0.4 any longer\n",
             "  because some of the arguments have changed")
    stopifnot(!missing(target))

    if (field < 0)
        stop("'field' must be an integer >=0")

    in.shapefile = default.file.extension(in.shapefile, ".shp")
    out.grid = default.file.extension(out.grid, ".sgrd")
    
    if (target$TARGET == 1) {
        if (target$GRID_GRID != out.grid) {
            rsaga.copy.sgrd(target$GRID_GRID, out.grid, env = env)
            target$GRID_GRID = out.grid
        }
    }

    param = list(
        USER_GRID = out.grid,
        SHAPES = in.shapefile,
        FIELD = field)
    param = c(param, target)
    
    # TARGET parameters changed SAGA 2.1.3:
    if (env$version == "2.1.3" | env$version == "2.1.4" | env$version == "2.2.0" |
        env$version == "2.2.1" | env$version == "2.2.2" | env$version == "2.2.3") {
        nm = names(param)
        nm[ nm == "USER_GRID" ] = "TARGET_OUT_GRID"
        nm[ nm == "TARGET" ] = "TARGET_DEFINITION"
        nm[ nm == "GRID_GRID" ] = "TARGET_TEMPLATE"
        nm[ nm == "USER_SIZE" ] = "TARGET_USER_SIZE"
        nm[ nm == "USER_FIT" ] = "TARGET_USER_FITS"
        nm[ nm == "USER_XMIN" ] = "TARGET_USER_XMIN"
        nm[ nm == "USER_XMAX" ] = "TARGET_USER_XMAX"
        nm[ nm == "USER_YMIN" ] = "TARGET_USER_YMIN"
        nm[ nm == "USER_YMAX" ] = "TARGET_USER_YMAX"
        names(param) = nm
    }

    rsaga.geoprocessor(lib = "grid_gridding", 
        module = "Nearest Neighbour", # was: = 2 (=1 in earlier SAGA version)
        param, env = env, ...)
}

#' @rdname rsaga.inverse.distance
#' @name rsaga.modified.quadratic.shephard
#' @export
rsaga.modified.quadratic.shephard = function(in.shapefile, out.grid, field,
    quadratic.neighbors = 13, weighting.neighbors = 19,
    target, env = rsaga.env(), ...)
{
    if (env$version == "2.0.4")
        stop("rsaga.modified.quadratic.shephard doesn't support SAGA GIS 2.0.4 any longer\n",
             "  because some of the arguments have changed")
    stopifnot(!missing(target))

    if (field < 0)
        stop("'field' must be an integer >=0")
    if (quadratic.neighbors < 5)
        stop("'quadratic.neighbors' must be an integer >=5")
    if (weighting.neighbors < 5)
        stop("'weighting.neighbors' must be an integer >=3")

    in.shapefile = default.file.extension(in.shapefile, ".shp")
    out.grid = default.file.extension(out.grid, ".sgrd")
    
    if (target$TARGET == 1) {
        if (target$GRID_GRID != out.grid) {
            rsaga.copy.sgrd(target$GRID_GRID, out.grid, env = env)
            target$GRID_GRID = out.grid
        }
    }

    param = list(
        USER_GRID = out.grid,
        SHAPES = in.shapefile,
        FIELD = field,
        QUADRATIC_NEIGHBORS = quadratic.neighbors,
        WEIGHTING_NEIGHBORS = weighting.neighbors)
    
    param = c(param, target)
    
    # TARGET parameters changed SAGA 2.1.3:
    if (env$version == "2.1.3" | env$version == "2.1.4" | env$version == "2.2.0" |
        env$version == "2.2.1" | env$version == "2.2.2" | env$version == "2.2.3") {
        nm = names(param)
        nm[ nm == "USER_GRID" ] = "TARGET_OUT_GRID"
        nm[ nm == "TARGET" ] = "TARGET_DEFINITION"
        nm[ nm == "GRID_GRID" ] = "TARGET_TEMPLATE"
        nm[ nm == "USER_SIZE" ] = "TARGET_USER_SIZE"
        nm[ nm == "USER_FIT" ] = "TARGET_USER_FITS"
        nm[ nm == "USER_XMIN" ] = "TARGET_USER_XMIN"
        nm[ nm == "USER_XMAX" ] = "TARGET_USER_XMAX"
        nm[ nm == "USER_YMIN" ] = "TARGET_USER_YMIN"
        nm[ nm == "USER_YMAX" ] = "TARGET_USER_YMAX"
        names(param) = nm
    }
        
    rsaga.geoprocessor(lib = "grid_gridding", 
        module = "Modifed Quadratic Shepard", # = 4 (earlier SAGA versions: =2)
        param, env = env, ...)
}


#' @rdname rsaga.inverse.distance
#' @name rsaga.triangulation
#' @export
rsaga.triangulation = function(in.shapefile, out.grid, field,
    target, env = rsaga.env(), ...)
{
    if (env$version == "2.0.4")
        stop("rsaga.triangulation doesn't support SAGA GIS 2.0.4 any longer\n",
             "  because some of the arguments have changed")
    stopifnot(!missing(target))

    if (field < 0)
        stop("'field' must be an integer >=0")

    in.shapefile = default.file.extension(in.shapefile, ".shp")
    out.grid = default.file.extension(out.grid, ".sgrd")
    
    if (target$TARGET == 1) {
        if (target$GRID_GRID != out.grid) {
            rsaga.copy.sgrd(target$GRID_GRID, out.grid, env = env)
            target$GRID_GRID = out.grid
        }
    }

    param = list(
        USER_GRID = out.grid,
        SHAPES = in.shapefile,
        FIELD = field)
    param = c(param, target)
    
    # TARGET parameters changed SAGA 2.1.3:
    if (env$version == "2.1.3" | env$version == "2.1.4" | env$version == "2.2.0" |
        env$version == "2.2.1" | env$version == "2.2.2" | env$version == "2.2.3") {
        nm = names(param)
        nm[ nm == "USER_GRID" ] = "TARGET_OUT_GRID"
        nm[ nm == "TARGET" ] = "TARGET_DEFINITION"
        nm[ nm == "GRID_GRID" ] = "TARGET_TEMPLATE"
        nm[ nm == "USER_SIZE" ] = "TARGET_USER_SIZE"
        nm[ nm == "USER_FIT" ] = "TARGET_USER_FITS"
        nm[ nm == "USER_XMIN" ] = "TARGET_USER_XMIN"
        nm[ nm == "USER_XMAX" ] = "TARGET_USER_XMAX"
        nm[ nm == "USER_YMIN" ] = "TARGET_USER_YMIN"
        nm[ nm == "USER_YMAX" ] = "TARGET_USER_YMAX"
        names(param) = nm
    }
        
    rsaga.geoprocessor(lib = "grid_gridding", 
        module = "Triangulation",
        param, env = env, ...)
}
