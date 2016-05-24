#' @title Run SaTScan in the OS
#' 
#' @description Calls out into the OS to run SaTScan, with the parameter file specified
#' 
#' @param prmlocation A string containing the 
#' directory location where the paramter file is located. 
#' @param prmfilename A string containg the name of the parameter file, without the extension, i.e., no ".prm".
#' @param sslocation A string containing the directory location where satscanbatch.exe (Windows) 
#' is located.  The default value is a common location in Windows 7.
#' @param ssbatchfilename Name of the file containing the SaTScan executable.  This is likely to be
#' either SaTScanBatch or SaTScanBatch64.  Omit the file extension.
#' @param cleanup If true, deletes any SaTScan output files from the OS.
#' @param verbose If true, will display the results in the R console as if running SaTScan in batch.  This may be especially
#' useful if you expect SaTScan to take a long time to run.
#' 
#' @details The parameter file may have been made by the \code{\link{ss.options}} function or not.  
#' If not, or if the \code{matchout = FALSE} parameter was set in \code{ss.options}, then the
#' return object will include the main text output from SaTScan only you manually set the 
#' \code{ResultsFile} SaTScan parameter to have the same name as the parameter file.
#' 
#' @return A satscan-class object, which is a list of 8 items, not all of which are always made, depending on SaTScan options and whether the program call was successful or not:
#' \describe{
#'   \item{main}{A character vector containing the main text output from SaTScan.  This is 
#'        probably identical to the material displayed when verbose=True}
#'   \item{col}{A data frame with the basic cluster information dataset SaTScan makes.}
#'   \item{rr}{A data frame with the risk ratio dataset SaTScan makes.}
#'   \item{gis}{A data frame with the geographic information dataset SaTScan makes.}
#'   \item{llr}{A data frame with the log likelihood ratios dataset SaTScan makes.}
#'   \item{sci}{A data frame with the other cluster information dataset SaTScan makes.}
#'   \item{shapeclust}{A list object, of class SpatialPolygonsDataFrame, defined by the \code{sp} 
#'   package.  It contains the ESRI shapefile(s) SaTScan makes.  This is made only if the \code{rgdal}
#'   package is available.}
#'   \item{prm}{A character vector containing the contents of the parameter file you told SaTScan 
#'   to use.}
#' }
#' If an item is not made by SaTScan, it will be NA.
#' 
#' @examples 
#' 
#' \dontrun{
#' ## Please see vignette("rsatscan"); example() code doesn't make sense since
#' ## all examples rely on calls to SaTScan in the OS.
#' }
#' 
#' 
#' @seealso \code{\link{ss.options}}, \code{\link{write.ss.prm}}
#' 
#' @export

satscan = function(
  prmlocation,  prmfilename,  sslocation = "c:/progra~2/satscan",
  ssbatchfilename = "SaTScanBatch", cleanup = TRUE, verbose=FALSE) {
    
    
    ssfile = paste0(stripslash(sslocation), "/", ssbatchfilename)
    if (Sys.which(ssfile) =="")
          stop("SaTScan is not there or is not runnable")
    prmloc = paste0(stripslash(prmlocation),"/")
    infile = paste0(prmloc, prmfilename,".prm")
    if (!file.exists(infile))  stop("I can't find that parameter file")
    system(paste(shQuote(ssfile), infile), show.output.on.console=verbose)
    prm = suppressWarnings(readLines(infile))
    mainfile = if  (file.exists(paste0(prmloc,prmfilename,".txt")))
      read.satscanmain(prmloc,prmfilename) else NA
    colfile = if  (file.exists(paste0(prmloc,prmfilename,".col.dbf")))
      read.col(prmloc,prmfilename) else NA
    rrfile = if  (file.exists(paste0(prmloc,prmfilename,".rr.dbf")))
      read.rr(prmloc,prmfilename) else NA
    gisfile = if  (file.exists(paste0(prmloc,prmfilename,".gis.dbf")))
      read.gis(prmloc,prmfilename) else NA
    llrfile = if  (file.exists(paste0(prmloc,prmfilename,".llr.dbf")))
      read.llr(prmloc,prmfilename) else NA
    scifile = if  (file.exists(paste0(prmloc,prmfilename,".sci.dbf")))
      read.sci(prmloc,prmfilename) else NA  
    if (verbose) cat("\n \n Any following message is from readOGR() in the rdgal package \n")
      shpfile = if (file.exists(paste0(prmloc,prmfilename,".col.shp"))) {
        if (requireNamespace("rgdal", quietly = TRUE)) {
          rgdal::readOGR(dsn=stripslash(prmloc), layer=paste0(prmfilename,".col"),
                    verbose=verbose) 
                }  else cat("\n \n rgdal package not installed, so shapefile can't be imported")
          } else NA
    if (cleanup) {
        suppressWarnings(file.remove(paste0(prmloc,prmfilename,".txt")))
        suppressWarnings(file.remove(paste0(prmloc,prmfilename,".col.dbf")))
        suppressWarnings(file.remove(paste0(prmloc,prmfilename,".rr.dbf")))
        suppressWarnings(file.remove(paste0(prmloc,prmfilename,".gis.dbf")))
        suppressWarnings(file.remove(paste0(prmloc,prmfilename,".llr.dbf")))
        suppressWarnings(file.remove(paste0(prmloc,prmfilename,".sci.dbf")))
        suppressWarnings(file.remove(paste0(prmloc,prmfilename,".col.shp")))
        suppressWarnings(file.remove(paste0(prmloc,prmfilename,".col.shx")))
      }
  return(structure(list(main=mainfile, col=colfile, rr = rrfile, gis=gisfile, llr=llrfile, 
                        sci=scifile, shapeclust=shpfile, prm=prm), class="satscan"))
}




