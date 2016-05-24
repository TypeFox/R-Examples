# Purpose        : Initial settings;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : Pierre Roudier (pierre.roudier@landcare.nz); Dylan Beaudette (debeaudette@ucdavis.edu); 
# Dev Status     : Pre-Alpha
# Note           : for more info see [http://cran.r-project.org/doc/manuals/R-exts.html]; this code was prepared for SAGA GIS 2.0.8


################## STANDARD ENVIRONMENTS ##############

## setup our environment for storing file handles and the like
plotKML.fileIO <- new.env(hash=TRUE, parent = parent.frame())

## setup the plotKML environment:
plotKML.opts <- new.env(hash=TRUE, parent = parent.frame())

## Find paths to external packages:
paths <- function(gdalwarp = "", gdal_translate = "", convert = "", saga_cmd = "", python = "", gdal.dir = utils::shortPathName("C:\\Program Files\\GDAL"), show.paths = TRUE){ 

     ## Try locating SAGA GIS (R default setting)...
     if(saga_cmd==""){
      #require(RSAGA)
      if(!inherits(try( suppressWarnings( x <- rsaga.env() ), silent = TRUE), "try-error")){
        if(!is.null(x)){ 
          if(.Platform$OS.type == "windows") {
            saga_cmd <- utils::shortPathName(normalizePath(paste(rsaga.env()$path, rsaga.env()$cmd, sep="/"))) 
          } else { 
            saga_cmd <- paste(rsaga.env()$path, rsaga.env()$cmd, sep="/") 
          } 
        if(nzchar(saga_cmd)){
          saga.version <- rsaga.get.version()
        }
      } else {
        saga.version <- ""
        }
      }
     }
     
     ## Try locating path to ImageMagick (R default setting)...
     if(convert==""){
       if(requireNamespace("animation", quietly = TRUE)){
          convert <- animation::ani.options("convert")
       }

     ## If it does not work, try getting the path from the OS:     
     if(is.null(convert)){
        im.dir <- NULL
        if(.Platform$OS.type == "windows") {
          ## get paths and check for ImageMagick
          paths <- strsplit(Sys.getenv('PATH')[[1]], ";")[[1]]
          x <- grep(paths, pattern="Magick", ignore.case = TRUE)
          
          ## if present
          if(!length(x) == 0) {
            im.dir <- paths[grep(paths, pattern="Magick", ignore.case = TRUE)[1]]
            convert = shQuote(utils::shortPathName(normalizePath(file.path(im.dir, "convert.exe"))))
            if(show.paths&file.exists(convert)){ 
              try( om <- system(convert,  show.output.on.console = FALSE, intern = TRUE)[1] )
              if(!class(.Last.value)[1]=="try-error"){
                message( paste(om) ) 
              } else {
                convert <- NULL
              }
            }
          }
        } ## end checking for Imagemagick on Windows
        
        ## check for all other OS:
        else {
          if(!length(x <- grep(paths <- strsplit(Sys.getenv('PATH')[[1]], ":")[[1]], pattern="Magick", ignore.case = TRUE))==0) {
            im.dir <- paths[grep(paths, pattern="Magick", ignore.case = TRUE)[1]]
            convert = "convert"
            if(show.paths&file.exists(convert)){ 
              try( om <- system(convert,  show.output.on.console = FALSE, intern = TRUE)[1] )
              if(!class(.Last.value)[1]=="try-error"){
                message( paste(om) ) 
              } else {
                convert <- NULL
              } 
            }
          } else {
            im.dir <- NULL
          }
        }
    
        if(is.null(im.dir)){ 
          warning("Install ImageMagick and add to PATH. See http://imagemagick.org for more info.")
        convert = ""
        }
     } else { 
     if(show.paths&file.exists(convert)){ 
       try( om <- system(convert,  show.output.on.console = FALSE, intern = TRUE)[1] )
       if(!class(.Last.value)[1]=="try-error"){
         message( paste(om) ) 
       } else {
         convert <- NULL
       } 
     }
    }
    }  
  
    ## try to locate GDAL / Patyhon:
    if(.Platform$OS.type == "windows") {
      if(gdalwarp==""|gdal_translate==""){
        if(requireNamespace("gdalUtils", quietly = TRUE)){
          gdalUtils::gdal_setInstallation(search_path=gdal.dir, rescan=FALSE)
          x <- getOption("gdalUtils_gdalPath")
          if(!is.null(x[[1]]$path)){
            gdalwarp = shQuote(utils::shortPathName(normalizePath(file.path(x[[1]]$path, "gdalwarp.exe"))))
            gdal_translate = shQuote(utils::shortPathName(normalizePath(file.path(x[[1]]$path, "gdal_translate.exe"))))
        }} else {
          warning("Could not locate GDAL! Install program and add it to the Windows registry. See http://www.gdal.org/ for more info.")
          gdalwarp = ""
          gdal_translate = ""    
        } 
      }
        
        ## Python:
        if(python==""){
        reg.paths <- names(utils::readRegistry("SOFTWARE"))
        x <- grep(reg.paths, pattern="WOW6432Node", ignore.case = TRUE)
        if(length(x)>0 & !inherits(try({ 
          py.paths <- utils::readRegistry(paste("SOFTWARE", reg.paths[x], "Python", sep="\\"), maxdepth=3)
          py.path = utils::readRegistry(paste("SOFTWARE", reg.paths[x], "Python", names(py.paths), names(py.paths[[1]]), "InstallPath", sep="\\"))[[1]] 
          }, silent = TRUE), "try-error")) {
            if (nzchar(py.path))  { 
              python = shQuote(utils::shortPathName(normalizePath(file.path(py.path, "python.exe"))))
              if(show.paths){ 
                message(paste("Located Python from the Registry Hive: \"", utils::shortPathName(py.path), "\"", sep="")) 
              }
            } 
        } else { 
        if(!inherits(try({ 
          py.paths <- utils::readRegistry(paste("SOFTWARE", "Python", sep="\\"), maxdepth=3)
          py.path = utils::readRegistry(paste("SOFTWARE", "Python", names(py.paths), names(py.paths[[1]])[1], "InstallPath", sep="\\"))[[1]] 
          }, silent = TRUE), "try-error")) {
          if (nzchar(py.path))  { 
            python = shQuote(utils::shortPathName(normalizePath(file.path(py.path, "python.exe"))))
            if(show.paths){ 
              message(paste("Located Python from the Registry Hive: \"", utils::shortPathName(py.path), "\"", sep="")) 
            }
          }
        } 
             
        }
        }
         
        ## 2nd chance to try to locate SAGA GIS (if not on a standard path):      
        if(saga_cmd==""){
          if(!nzchar(saga_cmd)&!nzchar(saga.version)){
          if(nzchar(prog <- Sys.getenv("ProgramFiles")) &&
            length(saga.dir <- list.files(prog, "^SAGA*"))>0 &&
            length(saga_cmd <- list.files(file.path(prog, saga.dir), pattern = "^saga_cmd\\.exe$", full.names = TRUE, recursive = TRUE))>0  
            ){
          if(suppressWarnings(!is.null(myenv <- rsaga.env(path=shQuote(normalizePath(saga.dir[1])))))){ 
            saga_cmd <- utils::shortPathName(normalizePath(paste(myenv$path, myenv$cmd, sep="/")))
            saga.version <- myenv$version 
            if(show.paths){ 
              message(paste("Located SAGA GIS ", saga.version, " from the 'Program Files' directory: \"", utils::shortPathName(saga_cmd), "\"", sep="")) 
            }
       }} else{ if(nzchar(prog <- Sys.getenv("ProgramFiles(x86)")) &&
            length(saga.dir <- list.files(prog, "^SAGA*"))>0 &&
            length(saga_cmd <- list.files(file.path(prog, saga.dir), pattern = "^saga_cmd\\.exe$", full.names = TRUE, recursive = TRUE))>0   
            ) {
          if(suppressWarnings(!is.null(myenv <- rsaga.env(path=shQuote(normalizePath(saga.dir[1])))))){ 
            saga_cmd <- utils::shortPathName(normalizePath(paste(myenv$path, myenv$cmd, sep="/")))
            saga.version <- myenv$version 
            if(show.paths){ 
              message(paste("Located SAGA GIS ", saga.version, " from the 'Program Files' directory: \"", utils::shortPathName(saga_cmd), "\"", sep="")) 
            }
       }}
       }
        
       if(!nzchar(saga_cmd)){
          warning("Could not locate SAGA GIS! Install program and add it to the Windows registry. See http://www.saga-gis.org/en/ for more info.")
          saga_vc = "" 
        }   
       }
       else {
          if(show.paths){ 
            message(paste("Located SAGA GIS ", saga.version, " from the 'Program Files' directory: \"", utils::shortPathName(saga_cmd), "\"", sep="")) 
          }
       }
      }
    }
    
    ## UNIX:
    else {
    
    if(gdalwarp==""|gdal_translate==""){
    if(!length(x <- grep(paths <- strsplit(Sys.getenv('PATH')[[1]], ":")[[1]], pattern="GDAL", ignore.case=TRUE))==0) {
    gdal.dir <- paths[grep(paths, pattern="GDAL", ignore.case=TRUE)[1]]
    gdalwarp = "gdalwarp"
    gdal_translate = "gdal_translate"
    if(show.paths){ 
      message(paste("Located GDAL from the path: \"", gdal.dir, "\"", sep=""))
    }
    }
    else { 
        warning("Install GDAL and add to PATH. See http://www.gdal.org/ for more info.")
      gdalwarp = ""
      gdal_translate = ""
    }
    }
    
    if(python==""){
    if(!length(x <- grep(paths <- strsplit(Sys.getenv('PATH')[[1]], ":")[[1]], pattern="Python", ignore.case=TRUE))==0) {
    py.dir <- paths[grep(paths, pattern="Python", ignore.case=TRUE)[1]] 
    python = "python"
    if(show.paths){ message(paste("Located Python from the path: \"", py.dir, "\"", sep="")) }
      }
    else { 
        warning("Install Python and add to PATH. See http://python.org for more info.")
        python = ""
    }
    }
    
    if(convert==""){
    im.dir <- paths[grep(paths, pattern="Magick", ignore.case=TRUE)[1]]
    if(is.null(im.dir)){ 
        warning("Install ImageMagick and add to PATH. See http://imagemagick.org for more info.")
        convert = ""
    }
    }

    if(saga_cmd==""){
    if(!nzchar(saga_cmd)){
        warning("Install SAGA GIS and add to PATH. See http://www.saga-gis.org for more info.")
        } 
    }
    }

    lt <- data.frame(gdalwarp, gdal_translate, convert, python, saga_cmd, stringsAsFactors = FALSE)
    return(lt)
}

################## STANDARD SETTINGS ##############

plotKML.env <- function(
    colour_scale_numeric = '',
    colour_scale_factor = '',
    colour_scale_svar = '',
    ref_CRS,
    NAflag,
    icon,
    LabelScale,
    size_range,
    license_url,
    metadata_sel,
    kmz,
    kml_xsd,
    kml_url,
    kml_gx,
    gpx_xsd,
    fgdc_xsd,
    inspire_xsd,
    convert,
    gdalwarp,
    gdal_translate,
    python,
    home_url,
    show.env = TRUE,
    silent = TRUE
    ){
    
	brewer1 = c("#D7191C","#FDAE61","#FFFFBF","#ABD9E9","#2C7BB6")
    #if(missing(colour_scale_numeric)) { colour_scale_numeric <- rev(brewer.pal(n = 5, name = "RdYlBu")) }
    if(missing(colour_scale_numeric)) { colour_scale_numeric <- rev(brewer1) }
	brewer2 = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999")
    #if(missing(colour_scale_factor)) { colour_scale_factor <- brewer.pal(n = 9, name = "Set1") }
    if(missing(colour_scale_factor)) { colour_scale_factor <- brewer2 }
	brewer3 = c("#FEEDDE","#FDBE85","#FD8D3C","#E6550D","#A63603")
    #if(missing(colour_scale_svar)) { colour_scale_svar <- brewer.pal(n = 5, name = "Oranges") }
    if(missing(colour_scale_svar)) { colour_scale_svar <- brewer3 }
    if(missing(ref_CRS)) { ref_CRS <- "+proj=longlat +datum=WGS84" }
    if(missing(NAflag)) { NAflag <- -99999 }
    if(missing(icon)) { icon <- "icon3.png" }   # "http://maps.google.com/mapfiles/kml/shapes/donut.png"
    if(missing(LabelScale)) { LabelScale <- .5 }
    if(missing(size_range)) { size_range <- c(0.25, 2.5) }
    if(missing(license_url)) { license_url <- "http://creativecommons.org/licenses/by/3.0/" }
    if(missing(metadata_sel)) { metadata_sel <- c("idinfo.citation.citeinfo.title", "idinfo.descript.abstract", "spdoinfo.ptvctinf.sdtsterm.ptvctcnt", "idinfo.timeperd.timeinfo.rngdates.begdate", "idinfo.timeperd.timeinfo.rngdates.enddate", "distinfo.stdorder.digform.digtopt.onlinopt.computer.networka.networkr", "idinfo.citation.citeinfo.othercit", "idinfo.citation.citeinfo.onlink", "idinfo.datacred", "distinfo.distrib.cntinfo.cntorgp.cntorg", "distinfo.stdorder.digform.digtinfo.formcont", "idinfo.native") }   
    if(missing(kmz)) { kmz <- FALSE }
    if(missing(kml_xsd)) { kml_xsd <- "http://schemas.opengis.net/kml/2.2.0/ogckml22.xsd" }
    if(missing(kml_url)) { kml_url <- "http://www.opengis.net/kml/2.2/" }
    if(missing(kml_gx)) { kml_gx <- "http://www.google.com/kml/ext/2.2" }
    if(missing(gpx_xsd)) { gpx_xsd <- "http://www.topografix.com/GPX/1/1/gpx.xsd" }
    if(missing(fgdc_xsd)) { fgdc_xsd <- "http://fgdcxml.sourceforge.net/schema/fgdc-std-012-2002/fgdc-std-012-2002.xsd" }
    if(missing(inspire_xsd)) { inspire_xsd <- "http://inspire.ec.europa.eu/schemas/common/1.0/common.xsd" }
    
    if(silent == FALSE){
      pts <- paths(show.paths = TRUE)
    }
    else {
      pts <- data.frame(gdalwarp="", gdal_translate="", convert="", python="", saga_cmd="", stringsAsFactors = FALSE)
    }
    
    if(missing(convert)) { convert <- pts$convert[[1]] }
    if(missing(gdalwarp)) { gdalwarp <- pts$gdalwarp[[1]] }
    if(missing(gdal_translate)) { gdal_translate <- pts$gdal_translate[[1]] }
    if(missing(python)) { python <- pts$python[[1]] }
    if(missing(home_url)) { home_url <- "http://plotkml.r-forge.r-project.org/" }
 
    ## Create a list and assing it to plotKML.env:
    pl.lst <- list(
      colour_scale_numeric = colour_scale_numeric,
      colour_scale_factor = colour_scale_factor,
      colour_scale_svar = colour_scale_svar,
      ref_CRS = ref_CRS,
      NAflag = NAflag,
      icon = icon,
      LabelScale = LabelScale,
      size_range = size_range,
      license_url = license_url,
      metadata_sel = metadata_sel,
      kmz = kmz,
      kml_xsd = kml_xsd,
      kml_url = kml_url,
      kml_gx = kml_gx,
      gpx_xsd = gpx_xsd,
      fgdc_xsd = fgdc_xsd,
      inspire_xsd = inspire_xsd,
      convert = convert,
      gdalwarp = gdalwarp,
      gdal_translate = gdal_translate,
      python = python,
      home_url = home_url
    )
 
    x <- lapply(names(pl.lst), function(x){ assign(x, pl.lst[[x]], envir=plotKML.opts) })
     
    if(show.env){  return(pl.lst)  }
 
}

## generate all environmental settings:
plotKML.env(show.env = FALSE)

# end of script;
