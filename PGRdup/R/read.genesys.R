### This file is part of 'PGRdup' package for R.

### Copyright (C) 2014, ICAR-NBPGR.
#
# PGRdup is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# PGRdup is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#' Convert 'Darwin Core - Germplasm' zip archive to a flat file
#' 
#' \code{read.genesys} reads PGR data in a Darwin Core - germplasm zip 
#' archive downloaded from genesys database and creates a flat file
#' \code{data.frame} from it.
#' 
#' This function helps to import to R environment, the PGR data 
#' downloaded from genesys database \url{https://www.genesys-pgr.org/} as a 
#' Darwin Core - germplasm (DwC-germplasm) zip archive. The different csv files 
#' in the archive are merged as a flat file into a single \code{data.frame}.
#' 
#' All the space characters can be removed from the fields corresponding to 
#' accession names such as acceNumb, collNumb, ACCENAME, COLLNUMB, DONORNUMB and
#' OTHERNUMB using the argument \code{scrub.names.space} to facilitate creation 
#' of KWIC index with \code{\link[PGRdup]{KWIC}} function and subsequent
#' matching operations to identify probable duplicates with
#' \code{\link[PGRdup]{ProbDup}} function.
#' 
#' The argument \code{readme} can be used to print the readme file in the 
#' archive to console, if required.
#' 
#' @param zip.genesys A character vector giving the file path to the downloaded 
#'   zip file from Genesys.
#' @param scrub.names.space logical. If \code{TRUE}, all space characters are 
#'   removed from name field in names extension (see \strong{Details}).
#' @param readme logical. If \code{TRUE}, the genesys zip file readme is printed
#'   to console.
#' @return A data.frame with the flat file form of the genesys data.
#' @seealso \code{\link[data.table]{data.table}}
#' 
#' @import data.table
#' @importFrom utils unzip
#' @export
read.genesys <- function(zip.genesys, scrub.names.space = TRUE, readme = TRUE) {
  # Check whether archive is a Darwin Core - Germplasm Archive
  FileList <- unzip(zip.genesys, list = TRUE)$Name
  DwCList <- c("README.txt", "core.csv", "names.csv", "geo.csv", "coll.csv",
               "meta.xml")
  if (setequal(FileList, DwCList) == FALSE) {
    stop("The zip file is not a Genesys Darwin Core - Germplasm Archive")
  }
  # Create the temporary directory or flush CSVs if it exists already
  if (!file.exists(tempdir())) {
    dir.create(tempdir())
  } else {
    files <- c(paste(tempdir(), "\\names.csv", sep = ""),
               paste(tempdir(), "\\geo.csv", sep = ""),
               paste(tempdir(), "\\core.csv", sep = ""),
               paste(tempdir(), "\\coll.csv", sep = ""))
    for (i in 1:4) {
      if (file.exists(files[i])) file.remove(files[i])
    }
  }
  # Unzip the file into the dir
  unzip(zip.genesys, exdir = tempdir())
  # Import names.csv
  nam <- fread(input = list.files(tempdir(), pattern = "names.csv",
                                  full.names = T),
               header = TRUE, stringsAsFactors = FALSE,
               select = c("genesysId", "instCode", "name", "aliasType",
                          "lang", "version"),
               colClasses = rep("character", 6))
  setkey(nam, genesysId)
  if ( !length(setdiff(nam$aliasType,
                      c("ACCENAME", "COLLNUMB",
                        "DONORNUMB", "OTHERNUMB"))) == 0) {
    warning("Abnormal strings detected in 'aliasType' column of 'names.csv'")
  }

  # Remove space from acc name fields 1
  if (scrub.names.space) {
    nam[, name := gsub(pattern = "[[:space:]]", replacement = "", name)]
  }
  # Convert to wide form
  nam <- dcast.data.table(nam, genesysId ~ aliasType, value.var = "name",
                          fun.aggregate = function(x) paste(x, collapse = ":"))
  if (!length(setdiff(colnames(nam),
                     c("genesysId", "ACCENAME", "COLLNUMB",
                       "DONORNUMB", "OTHERNUMB"))) == 0) {
    nam[, setdiff(colnames(nam),
                  c("genesysId", "ACCENAME", "COLLNUMB",
                    "DONORNUMB", "OTHERNUMB")) := NULL]
  }

  # Import geo.csv
  geo <- fread(input = list.files(tempdir(), pattern = "geo.csv",
                                  full.names = T),
               header = TRUE, stringsAsFactors = FALSE,
               select = c("genesysId", "latitude", "longitude", "elevation",
                          "datum", "uncertainty", "method", "version"),
               colClasses = rep("character", 8))
  setkey(geo, genesysId)
  # Import core.csv
  cor <- fread(input = list.files(tempdir(), pattern = "core.csv",
                                  full.names = T),
               header = TRUE, stringsAsFactors = FALSE,
               select = c("genesysId", "uuid", "instCode", "acceNumb", "genus",
                          "species", "fullTaxa", "orgCty", "acqSrc", "acqDate",
                          "mlsStat", "available", "historic", "storage",
                          "sampStat", "duplSite", "createdBy", "createdDate",
                          "lastModifiedBy", "lastModifiedDate"),
               colClasses = rep("character", 20))
  setkey(cor, genesysId)
  # Import coll.csv
  col <- fread(input = list.files(tempdir(), pattern = "coll.csv",
                                  full.names = T),
               header = TRUE, stringsAsFactors = FALSE,
               select = c("genesysId", "collMissId", "collNumb", "collDate",
                          "collSrc", "collSite", "collCode", "collName",
                          "collInstAddress", "version"),
               colClasses = structure(rep("character", 10),
                                      .Names = c("genesysId", "collMissId",
                                                 "collNumb", "collDate",
                                                 "collSrc", "collSite",
                                                 "collCode", "collName",
                                                 "collInstAddress", "version")))
  setkey(col, genesysId)
  # Import and print the readme
  con <- unz(zip.genesys, "README.txt")
  rdm <- readLines(con)
  on.exit(close(con))
  if (readme) {
    cat(rdm, sep = "\n")
  }
  # Merge all csv files
  m <- merge(cor, col, by = "genesysId", all = TRUE)
  m <- merge(m, nam, by = "genesysId", all = TRUE)
  if (dim(geo)[1] != 0) {
    m <- merge(m, geo, by = "genesysId", all = TRUE,
               suffixes = c("_COLL", "_GEO"))
  }
  # Remove space from acc name fields 2
  if (scrub.names.space) {
    ts <-  c("acceNumb", "collNumb")
    m[, (ts) := lapply(.SD, function(x) gsub(pattern = "[[:space:]]",
                       replacement = "", x)), .SDcols = ts]
  }
  m <- as.data.frame(m)
  attr(m, "readme") <- rdm
  return(m)
  }
