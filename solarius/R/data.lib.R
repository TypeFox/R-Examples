#----------------------------------
# Data loaders
#----------------------------------

#' Alternative to system.file
#'
#' The function works as \code{system.file},
#' but takes care when the package is a local folder.
#' 
#' The use case is when some data file or direcotry is needed to be loaded,
#' and it is placed in \code{inst/} directory of the package.
#'
#' @param ...
#'    arguments to be passed to \code{system.file}
#' @return
#'    Path returned by \code{system.file}.
#'
#' @examples
#' mibddir <- package.file("extdata", "solarOutput", "solarMibds", package = "solarius") 
#' mibddir
#'
#' list.files(mibddir)
#'
#' @export
package.file <- function(...)
{
  path <- system.file(...)
  if(!file.exists(path)) {
    path <- system.file("inst", ...)
  }
  stopifnot(file.exists(path))
  
  return(path)
}

#' Load the complete data set from R package multic
#'
#' The function loads the complete data of 12,000 individuals,
#' which is stored in .phen and .ped files.
#' These two files were generated within R package \code{multic}
#' and re-distributed in R package \code{solarius} 
#' (\code{extdata/solarOutput} directory).
#'
#' Function \code{\link{readPhen}} is used to read .phen and .ped files.
#'
#' @examples
#' dat <- loadMulticPhen()
#' dim(dat)
#' 
#' data(dat30)
#' dim(dat30)
#' @export
loadMulticPhen <- function()
{
  dat.dir <- package.file("extdata", "solarOutput", package = "solarius")
  
  dat <- readPhen(phen.file = file.path(dat.dir, "simulated.phen"), sep.phen = ",", 
    ped.file = file.path(dat.dir, "simulated.ped"), sep.ped = ",")
  
  return(dat)
}


#' Load the example data (GAW10) from SOLAR.
#'
#' Function loads the example data distributed with SOLAR.
#' See the SOLAR tutorial for description of the data set.
#'
#' @param dat.dir
#'  (optional) A character with path to a directory, where SOLAR files to be exported.
#'  If this argument is missing, a temporary directory is used.
#' @return
#'  A data frame returned by \code{readPhen} function called
#'  with phen and ped files of the data set.
#'
#' @export
loadExamplesPhen <- function(dat.dir)
{ 
  ### dir
  is.tmpdir <- missing(dat.dir)
  
  if(is.tmpdir) {
    dir <- tempfile(pattern = "loadExamplesPhen-")
    stopifnot(dir.create(dir))
  }
  
  ### run solar
  if(is.tmpdir) {
    ret <- solar("example", dir)
    dat.dir <- dir
  }
  
  ### read files
  dat <- readPhen(phen.file = file.path(dat.dir, "gaw10.phen"), sep.phen = ",", 
    ped.file = file.path(dat.dir, "gaw10.ped"), sep.ped = ",")
  
  ### clean
  if(is.tmpdir) {
    unlink(dir, recursive = TRUE)
  }
  
  return(dat)
}


