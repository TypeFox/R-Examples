.onLoad <- function(libname, pkgname) {
  # Wickham it's the man!
  op <- options()
  op.water <- list(
    waterOverwrite = TRUE,
    waterWriteResults = FALSE,
    waterOutputFolder = ".",
    waterSRTMrepo = NULL,
    waterAutoAoi = TRUE
  )
  toset <- !(names(op.water) %in% names(op))
  if(any(toset)) options(op.water[toset])
  invisible()
}

# .onAttach <- function(libname, pkgname) {
#   packageStartupMessage("This package writes function's results to working directory. You can change")
#   packageStartupMessage("output folder, or completely disable this feature using waterOptions()")
# }


#' Global options for water package
#' @description 
#' This function is based on raster::rasterOptions by Robert Hijmans. 
#' @param overwrite     Logical. If TRUE and writeResults is TRUE it will 
#' overwrite results. If FALSE, results are save with a name with name_datetime.
#' @param writeResults  Logical. If TRUE it'll write result to disk. This is 
#' slower but if FALSE you can have out-of-memory problems.
#' @param outputFolder  Name of a folder to save files, relative to workind
#' folder. 
#' @param SRTMrepo      A folder where SRTM grids are stored, to create DEM. See
#' prepareSRTMdata()
#' @param autoAoi       Logical. If TRUE it'll look for a object called aoi on 
#' .GlobalEnv and use it as aoi. See createAoi()
#' @param default       Logical. If TRUE will revert all options to defaults 
#' values
#' @return 
#' list of the current options (invisibly). If no arguments are provided the options are printed.
#' @references 
#' Robert J. Hijmans (2015). raster: Geographic Data Analysis and Modeling. R
#' package version 2.4-18. http://CRAN.R-project.org/package=raster
#' @author Guillermo F Olmedo, \email{guillermo.olmedo@@gmail.com}
#' @export
waterOptions <- function (overwrite, writeResults, outputFolder,
                          SRTMrepo, autoAoi, default = FALSE) 
{
#   setFiletype <- function(format) {
#     if (.isSupportedFormat(format)) {
#       options(rasterFiletype = format)
#     }
#     else {
#       warning(paste("Cannot set filetype to unknown or unsupported file format:", 
#                     format, ". See writeFormats()"))
#     }
#   }
  setOverwrite <- function(overwrite) {
    if (is.logical(overwrite)) {
      options(waterOverwrite = overwrite)
    }
    else {
      warning(paste("Could not set overwrite. It must be a logical value"))
    }
  }
  setWriteResults <- function(write) {
    if (is.logical(write)) {
      options(waterWriteResults = write)
    }
    else {
      warning(paste("Could not set writeResults. It must be a logical value"))
    }
  }
  setOutputFolder <- function(tmpdir) {
    if (!missing(tmpdir)) {
      tmpdir <- trim(tmpdir)
      if (tmpdir != "") {
        lastchar = substr(tmpdir, nchar(tmpdir), nchar(tmpdir))
        if (lastchar != "/" & lastchar != "\\") {
          tmpdir <- paste0(tmpdir, "/")
        }
        options(waterOutputFolder = tmpdir)
      }
    }
  }
  setSRTMrepo <- function(srtmDir) {
    if (!missing(srtmDir)) {
      srtmDir <- trim(srtmDir)
      if (srtmDir != "") {
        lastchar = substr(srtmDir, nchar(srtmDir), nchar(srtmDir))
        if (lastchar != "/" & lastchar != "\\") {
          srtmDir <- paste0(srtmDir, "/")
        }
        options(waterSRTMrepo = srtmDir)
      }
    }
  }
  setAutoAoi <- function(autoAoi) {
    if (is.logical(autoAoi)) {
      options(waterAutoAoi = autoAoi)
    }
    else {
      warning(paste("Could not set autoAoi. It must be a logical value"))
    }
  }
#   setToDisk <- function(todisk) {
#     if (is.logical(todisk)) {
#       options(rasterToDisk = todisk)
#     }
#     else {
#       warning(paste("todisk argument must be a logical value"))
#     }
#   }
#   depracatedWarnings <- function(x) {
#     if (is.logical(x)) {
#       if (is.na(x)) {
#         x <- TRUE
#       }
#       options(rasterDepracatedWarnings = x)
#     }
#   }
  cnt <- 0
  if (default) {
    cnt <- 1
    options(waterOverwrite = TRUE)
    options(waterWriteResults = FALSE)
    options(waterOutputFolder = ".")
    options(waterSRTMrepo = NULL)
    options(waterAutoAoi = TRUE)
    v <- utils::packageDescription("water")[["Version"]]
  }
  if (!missing(overwrite)) {
    setOverwrite(overwrite)
    cnt <- cnt + 1
  } else {overwrite <- getOption("waterOverwrite")}
  if (!missing(writeResults)) {
    setWriteResults(writeResults)
    cnt <- cnt + 1
  } else {writeResults <- getOption("waterWriteResults")}
  if (!missing(outputFolder)) {
    setOutputFolder(outputFolder)
    cnt <- cnt + 1
  } else {outputFolder <- getOption("waterOutputFolder")}
  if (!missing(SRTMrepo)) {
    setSRTMrepo(SRTMrepo)
    cnt <- cnt + 1
  } else {SRTMrepo <- getOption("waterSRTMrepo")}
  if (!missing(autoAoi)) {
    setAutoAoi(autoAoi)
    cnt <- cnt + 1
  } else {autoAoi <- getOption("waterAutoAoi")}
  lst <- list(overwrite = overwrite, writeResults = writeResults,
              outputFolder = outputFolder, SRTMrepo= SRTMrepo,
              autoAoi = autoAoi)
  save <- FALSE
  if (save) {
    v <- utils::packageDescription("water")[["Version"]]
    fn <- paste(options("startup.working.directory"), "/rasterOptions_", 
                v, sep = "")
    oplst <- NULL
    oplst <- c(oplst, paste("rasterFiletype='", lst$format, 
                            "'", sep = ""))
    oplst <- c(oplst, paste("rasterOverwrite=", lst$overwrite, 
                            sep = ""))
    oplst <- c(oplst, paste("rasterDatatype='", lst$datatype, 
                            "'", sep = ""))
    oplst <- c(oplst, paste("rasterTmpDir='", lst$tmpdir, 
                            "'", sep = ""))
    oplst <- c(oplst, paste("rasterTmpTime='", lst$tmptime, 
                            "'", sep = ""))
    r <- try(write(unlist(oplst), fn), silent = TRUE)
    cnt <- 1
  }
  if (cnt == 0) {
    cat("overwrite     :", lst$overwrite, "\n")
    cat("writeResults  :", lst$writeResults, "\n")
    cat("outputFolder  :", lst$outputFolder, "\n")
    cat("SRTMrepo      :", lst$SRTMrepo, "\n")
    cat("autoAoi       :", lst$autoAoi, "\n")
  }
  invisible(lst)
}


